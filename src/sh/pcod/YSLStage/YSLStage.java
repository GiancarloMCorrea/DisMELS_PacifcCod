/**
 * YSLStage.java
 *
 * Updates:
 * 20181011  1. Corrected gL formula.
 *           2. Removed fcnDevelopment to match Parameters class.
 *           3. Added "attached" as new attribute (necessary with updated DisMELS).
 *           4. Removed "diam" since it's replaced by "length"
 * 20190722: 1. Removed fields associated with egg stage attributes "devStage" and "density"
 * 20210209: 1. Converted to using IBMFunctions, added STDG functions, renamed 'length' to std_len.
 *           2. Added dry_wgt, converted a number of other variables to attributes
 *           3. Changed criteria for PNR (now similar to YSA).
 * 20210209: 1. Made changes to setInfoFromIndividual, setInfoFromSuperIndividual, setAttributes, 
 *              and clone methods to deal with differences in attributes coming from EggStage instances.
 *
 */

package sh.pcod.YSLStage;

import com.vividsolutions.jts.geom.Coordinate;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.List;
import org.openide.util.lookup.ServiceProvider;
import sh.pcod.EggStage.EggStageAttributes;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateSL;
import wts.models.DisMELS.IBMFunctions.Mortality.ConstantMortalityRate;
import wts.models.DisMELS.IBMFunctions.Mortality.InversePowerLawMortalityRate;
import wts.models.DisMELS.IBMFunctions.Movement.DielVerticalMigration_FixedDepthRanges;
import wts.models.DisMELS.IBMFunctions.SwimmingBehavior.ConstantMovementRateFunction;
import wts.models.DisMELS.framework.*;
import wts.models.DisMELS.framework.IBMFunctions.IBMFunctionInterface;
import wts.models.utilities.CalendarIF;
import wts.models.utilities.DateTimeFunctions;
import wts.roms.model.Interpolator3D;
import wts.roms.model.LagrangianParticle;

/**
 * Life stage class for Pacific cod yolk-sac larvae (YSL).
 * 
 * @author William Stockhausen
 * @author Sarah Hinckley
 */
@ServiceProvider(service=LifeStageInterface.class)
public class YSLStage extends AbstractLHS {
    
    /** flag to use Sarah's approach to first feeding */
    public static boolean useFirstFeedingSH = true;//TODO: should be a parameter?
    
        //Static fields    
            //  Static fields new to this class
    /* flag to do debug operations */
    public static boolean debug = false;
    /* Class for attributes SH_NEW */
    public static final String attributesClass = 
            sh.pcod.YSLStage.YSLStageAttributes.class.getName();
    /* Class for parameters */
    public static final String parametersClass = 
            sh.pcod.YSLStage.YSLStageParameters.class.getName();
    /* Class for feature type for point positions */
    public static final String pointFTClass = 
            wts.models.DisMELS.framework.LHSPointFeatureType.class.getName();
    /* Classes for next LHS SH_NEW Add fdl*/
    public static final String[] nextLHSClasses = new String[]{ 
            sh.pcod.YSLStage.YSLStage.class.getName(),
            sh.pcod.FDLStage.FDLStage.class.getName()};

    /* Classes for spawned LHS */
    public static final String[] spawnedLHSClasses = new String[]{};
    
    /* string identifying environmental field with copepod densities */
    private static final String Cop = "Cop";
    /* string identifying environmental field with euphausiid densities */
    private static final String Eup = "Eup";
    /* string identifying environmental field with neocalanus densities */
    private static final String NCa = "NCa";
    
    //Instance fields
            //  Fields hiding ones from superclass
    /* life stage atrbutes object */
    protected YSLStageAttributes atts = null;
    /* life stage parameters object */
    protected YSLStageParameters params = null;
    
    //  Fields new to class
        //fields that reflect parameter values
    protected boolean isSuperIndividual;
    protected double  horizRWP;
    protected double  minStageDuration;
    protected double  maxStageDuration;
    protected double  minStageSize;
    protected double  stageTransRate;
    protected boolean useRandomTransitions;
    
        //fields that reflect (new) attribute values
    /** flag indicating individual is attached to bottom */
    protected boolean attached = false;
    /** standard length (mm) */
    protected double std_len = 0;
    /** dry weight (mg) */
    protected double dry_wgt = 0;
    /** growth rate for standard length (mm/d) */
    protected double grSL = 0;
    /** growth rate for dry weight (1/d) */
    protected double grDW = 0;
    /** in situ temperature (deg C) */
    protected double temperature = 0;
    /** in situ salinity */
    protected double salinity = 0;
    /** in situ water density */
    protected double rho = 0;
   /** in situ copepod density (mg/m^3, dry wt) */
    protected double copepods = 0;
     /** in situ euphausiid density (mg/m^3, dry wt) */
    protected double euphausiids = 0;
     /** in situ neocalanoid density (mg/m^3, dry wt) */
    protected double neocalanus = 0;
    /** integrated criterion for yolk-sac absorption */
    protected double progYSA = 0;
    /** integrated criterion for point-of-no return */
    protected double progPNR = 0;
    
            //other fields
    /** number of individuals transitioning to next stage */
    private double numTrans;  
    protected double durPNR; //time (days) to point-of-no return based on current temperature
    protected double durYSA; //time (days) to yolk sac absorption based on current temperature
    protected double ageYSA; //age at which yolk-sac absorption occurred
    
    //initialized values
    protected final double rndFeed = Math.random(); //random value of cumulative probability at which feeding occurs
    protected double  prFeed       = 0.0;   //cumulative probability with time of first feeding (Sarah's approach)
    protected double  prNotFed     = 1.0;   //cumulative probability of NOT having fed
    protected double  indivCopWgt  = 1.0e-6;//typical weight for individual small copepod (kg)
    protected double  fCumHazFcn   = 0.0;   //cumulative hazard function for first feeding
    protected boolean hasFed       = false; //feeding flag
    
    /** IBM function selected for mortality */
    private IBMFunctionInterface fcnMortality = null; 
    /** IBM function selected for growth in SL */
    private IBMFunctionInterface fcnGrSL = null; 
    /** IBM function selected for growth in DW */
    private IBMFunctionInterface fcnGrDW = null; 
    /** IBM function selected for vertical movement */
    private IBMFunctionInterface fcnVM = null; 
    /** IBM function selected for vertical velocity */
    private IBMFunctionInterface fcnVV = null; 
    /** IBM function selected for time to point-of-no-return */
    private IBMFunctionInterface fcnPNR = null; 
    /** IBM function selected for time to yolk-sac absorption */
    private IBMFunctionInterface fcnYSA = null; 
    
    private int typeMort = 0;//integer indicating mortality function
    private int typeGrSL = 0;//integer indicating SL growth function
    private int typeGrDW = 0;//integer indicating DW growth function
    private int typeVM   = 0;//integer indicating vertical movement function
    private int typeVV   = 0;//integer indicating vertical velocity function
    private int typePNR  = 0;//integer indicating PNR function
    private int typeYSA  = 0;//integer indicating YSA function
    
    private static final Logger logger = Logger.getLogger(YSLStage.class.getName());
    
    /**
     * Creates a new instance of YSLStage.  
     *  This constructor should be used ONLY to obtain
     *  the class names of the associated classes.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     */
    public YSLStage() {
        super("");
        super.atts = atts;
        super.params = params;
    }
    
    /**
     * Creates a new instance of YSLStage with the given typeName.
     * A new id number is calculated in the superclass and assigned to
     * the new instance's id, parentID, and origID. 
     * 
     * The attributes are vanilla.  Initial attribute values should be set,
     * then initialize() should be called to initialize all instance variables.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     */
    public YSLStage(String typeName) 
                throws InstantiationException, IllegalAccessException {
        super(typeName);
        atts = new YSLStageAttributes(typeName);
        atts.setValue(YSLStageAttributes.PROP_id,id);
        atts.setValue(YSLStageAttributes.PROP_parentID,id);
        atts.setValue(YSLStageAttributes.PROP_origID,id);
        params = (YSLStageParameters) LHS_Factory.createParameters(typeName);
        super.atts = atts;
        super.params = params;
        setParameters(params);
    }

    /**
     * Creates a new instance of LHS with type name and
     * attribute values given by input String array.
     * 
     * Side effects:
     *  1. Calls createInstance(LifeStageAttributesInterface), with associated effects,
     *  based on creating an attributes instance from the string array.
     * /
     * @param strv - attributes as string array
     * @return - instance of LHS
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     */
    @Override
    public YSLStage createInstance(String[] strv) 
                        throws InstantiationException, IllegalAccessException {
        LifeStageAttributesInterface theAtts = LHS_Factory.createAttributes(strv);
        YSLStage lhs = createInstance(theAtts);
        return lhs;
    }

    /**
     * Creates a new instance of this LHS with attributes (including type name) 
     * corresponding to the input attributes instance.
     * 
     * Side effects:
     *  1. If theAtts id attribute is "-1", then a new (unique) id value is created 
     *  for the new LHS instance.
     *  2. If theAtts parentID attribute is "-1", then it is set to the value for id.
     *  3. If theAtts origID attribute is "-1", then it is set to the value for id.
     *  4. initialize() is called to initialize variables and convert position
     *   attributes.
     * /
     * @param theAtts - attributes instance
     * @return - instance of LHS
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     */
    @Override
    public YSLStage createInstance(LifeStageAttributesInterface theAtts)
                        throws InstantiationException, IllegalAccessException {
        YSLStage lhs = null;
        if (theAtts instanceof YSLStageAttributes) {
            lhs = new YSLStage(theAtts.getTypeName());
            long newID = lhs.id;//save id of new instance
            lhs.setAttributes(theAtts);
            if (lhs.atts.getID()==-1) {
                //constructing new individual, so reset id values to those of new
                //SH_NEW
                lhs.id = newID;
                lhs.atts.setValue(YSLStageAttributes.PROP_id,newID);
            }
            newID = (Long) lhs.atts.getValue(YSLStageAttributes.PROP_parentID);
            if (newID==-1) {
                lhs.atts.setValue(YSLStageAttributes.PROP_parentID,newID);
            }
            newID = (Long) lhs.atts.getValue(YSLStageAttributes.PROP_origID);
            if (newID==-1) {
                lhs.atts.setValue(YSLStageAttributes.PROP_origID,newID);
            }
        }
        lhs.initialize();//initialize instance variables
        return lhs;
    }

    /**
     *  Returns the associated attributes.  
     * @return associated YSLStageAttributes instance
     */
    @Override
    public YSLStageAttributes getAttributes() {
        return atts;
    }

    /**
     * Sets the values of the associated attributes object to those in the input
     * String[]. This does NOT change the typeNameof the LHS instance (or the 
     * associated LHSAttributes instance) on which the method is called.
     * Attribute values are set using SimpleBenthicLHSAttributes.setValues(String[]).
     * Side effects:
     *  1. If th new id attribute is not "-1", then its value for id replaces the 
     *      current value for the lhs.
     *  2. If the new parentID attribute is "-1", then it is set to the value for id.
     *  3. If the new origID attribute is "-1", then it is set to the value for id.
     *  4. initialize() is called to initialize variables and convert position
     *   attributes.
     * /
     * @param strv - attribute values as String[]
     */
    @Override
    public void setAttributes(String[] strv) {
        long aid;
        atts.setValues(strv);
        aid = atts.getValue(YSLStageAttributes.PROP_id, id);
        if (aid==-1) {
            //change atts id to lhs id
            atts.setValue(YSLStageAttributes.PROP_id, id);
        } else {
            //change lhs id to atts id
            id = aid;
        }
        aid = atts.getValue(YSLStageAttributes.PROP_parentID, id);
        if (aid==-1) {
            atts.setValue(YSLStageAttributes.PROP_parentID, id);
        }
        aid = atts.getValue(YSLStageAttributes.PROP_origID, id);
        if (aid==-1) {
            atts.setValue(YSLStageAttributes.PROP_origID, id);
        }
        initialize();//initialize instance variables
    }

    /**
     * Sets the attributes for the instance by copying values from the input.
     * This does NOT change the typeName of the LHS instance (or the associated 
     * LHSAttributes instance) on which the method is called.
     * 
     * Note that id, parentID, and origID are copied, as are other attributes equivalent
     * in both life stages. Attributes in this stage that are not present in newAtts
     * are added appropriately.
     * 
     *  Side effects:
     *      updateVariables() is called to update instance variables.
     *      Instance field "id" is also updated.
     * @param newAtts - should be instance of YSLStageAttributes or EggStageAttributes
     */
    @Override
    public void setAttributes(LifeStageAttributesInterface newAtts) {
        if (newAtts instanceof YSLStageAttributes) {
            YSLStageAttributes oldAtts = (YSLStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));
        } else if (newAtts instanceof EggStageAttributes) {
            EggStageAttributes oldAtts = (EggStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));
            //need to map attributes with different names correctly
            atts.setValue(YSLStageAttributes.PROP_SL,oldAtts.getValue(EggStageAttributes.PROP_SL, std_len));
            atts.setValue(YSLStageAttributes.PROP_DW,oldAtts.getValue(EggStageAttributes.PROP_DW, dry_wgt));
            atts.setValue(YSLStageAttributes.PROP_grSL,oldAtts.getValue(EggStageAttributes.PROP_grSL, grSL));
            atts.setValue(YSLStageAttributes.PROP_grDW,oldAtts.getValue(EggStageAttributes.PROP_grDW, grDW));
            //need to set attributes NOT included in EggStageAttributes
            //set prey concentrations based on current location
            double[] pos = lp.getIJK();
            copepods    = i3d.interpolateValue(pos,Cop,Interpolator3D.INTERP_VAL);
            euphausiids = i3d.interpolateValue(pos,Eup,Interpolator3D.INTERP_VAL);
            neocalanus  = i3d.interpolateValue(pos,NCa,Interpolator3D.INTERP_VAL);
            atts.setValue(YSLStageAttributes.PROP_copepod,copepods);
            atts.setValue(YSLStageAttributes.PROP_euphausiid,euphausiids);
            atts.setValue(YSLStageAttributes.PROP_neocalanus,neocalanus);
            //set the following attributes to initial values 
            atts.setValue(YSLStageAttributes.PROP_progYSA,progYSA);
            atts.setValue(YSLStageAttributes.PROP_progPNR,progPNR);
            atts.setValue(YSLStageAttributes.PROP_prNotFed,prNotFed);            
        } else {
            //TODO: should throw an error here
            logger.info("setAttributes(): no match for attributes type:"+newAtts.toString());
        }
        id = atts.getValue(YSLStageAttributes.PROP_id, id);
        //SH_NEW
        //updateVariables(); //TODO: check if this is necessary--it's also called in setInfoFromIndividual
    }
    
    /**
     *  Sets the associated attributes object. Use this after creating an LHS instance
     * as an "output" from another LHS that is functioning as an ordinary individual.
     * @param oldLHS
     */
    @Override
    public void setInfoFromIndividual(LifeStageInterface oldLHS){
        /** 
         * Since this is a single individual making a transition, we need to:
         *  1) copy the Lagrangian Particle from the old LHS
         *  2) start a new track from the current position for the oldLHS
         *  3) copy the attributes from the old LHS (id's should remain as for old LHS)
         *  4) set values for attributes NOT included in oldLHS
         *  5) set age in stage = 0
         *  6) set active and alive to true
         *  7) update local variables
         */
        //copy LagrangianParticle information
        this.setLagrangianParticle(oldLHS.getLagrangianParticle());
        //start track at last position of oldLHS track
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_PROJECTED),COORDINATE_TYPE_PROJECTED);
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_GEOGRAPHIC),COORDINATE_TYPE_GEOGRAPHIC);
        
        LifeStageAttributesInterface oldAtts = oldLHS.getAttributes();            
        setAttributes(oldAtts);
        
        //reset some attributes
        atts.setValue(YSLStageAttributes.PROP_ageInStage, 0.0);//reset age in stage
        atts.setValue(YSLStageAttributes.PROP_active,true);    //set active to true
        atts.setValue(YSLStageAttributes.PROP_alive,true);     //set alive to true

        id = atts.getID(); //reset id for current LHS to one from old LHS

        //update local variables to capture changes made here
        updateVariables();
    }
    
    /**
     *  Sets the associated attributes object. Use this after creating an LHS instance
     * as an "output" from another LHS that is functioning as a super individual.
     * @param oldLHS
     * @param numTrans
     */
    @Override
    public void setInfoFromSuperIndividual(LifeStageInterface oldLHS, double numTrans) {
        /** 
         * Since the old LHS instance is a super individual, only a part 
         * (numTrans) of it transitioned to the current LHS. Thus, we need to:
         *          1) copy the Lagrangian Particle from the old LHS
         *          2) start a new track from the current position for the oldLHS
         *          3) copy some attribute values from old stage and determine values for new attributes
         *          4) make sure id for this LHS is retained, not changed
         *          5) assign old LHS id to this LHS as parentID
         *          6) copy old LHS origID to this LHS origID
         *          7) set number in this LHS to numTrans
         *          8) reset age in stage to 0
         *          9) set active and alive to true
         *         10) update local variables to match attributes
         */
        //copy LagrangianParticle information
        this.setLagrangianParticle(oldLHS.getLagrangianParticle());
        //start track at last position of oldLHS track
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_PROJECTED),COORDINATE_TYPE_PROJECTED);
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_GEOGRAPHIC),COORDINATE_TYPE_GEOGRAPHIC);
        
        //copy some variables that should not change
        long idc = id;
        
        //copy the old attribute values
        LifeStageAttributesInterface oldAtts = oldLHS.getAttributes();            
        setAttributes(oldAtts);
        
        //reset some attributes and variables
        id = idc;
        atts.setValue(YSLStageAttributes.PROP_id,idc);//reset id to one for current LHS
        atts.setValue(YSLStageAttributes.PROP_parentID,
                      oldAtts.getValue(LifeStageAttributesInterface.PROP_id));//copy old id to parentID
        atts.setValue(YSLStageAttributes.PROP_number, numTrans);//set number to numTrans
        atts.setValue(YSLStageAttributes.PROP_ageInStage, 0.0); //reset age in stage
        atts.setValue(YSLStageAttributes.PROP_active,true);     //set active to true
        atts.setValue(YSLStageAttributes.PROP_alive,true);      //set alive to true
            
        //update local variables to capture changes made here
        updateVariables();
    }

    /**
     *  Returns the associated parameters.  
     * @return the associated YSLStageParameters instance
     */
    @Override
    public YSLStageParameters getParameters() {
        return params;
    }

    /**
     * Sets the parameters for the instance to a cloned version of the input.
     * @param newParams - should be instance of YSLStageParameters
     */
    @Override
    public void setParameters(LifeStageParametersInterface newParams) {
        if (newParams instanceof YSLStageParameters) {
            params = (YSLStageParameters) newParams;
            super.params = params;
            setParameterValues();
            fcnMortality = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_Mortality);
            fcnGrSL = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_GrowthSL);
            fcnGrDW = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_GrowthDW);
            fcnVM   = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_VerticalMovement);
            fcnVV   = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_VerticalVelocity);
            fcnPNR  = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_PNR);
            fcnYSA  = params.getSelectedIBMFunctionForCategory(YSLStageParameters.FCAT_YSA);
            
            if (fcnMortality instanceof ConstantMortalityRate)
                typeMort = YSLStageParameters.FCN_Mortality_ConstantMortalityRate;
            else if (fcnMortality instanceof InversePowerLawMortalityRate)
                typeMort = YSLStageParameters.FCN_Mortality_InversePowerLawMortalityRate;
            
            if (fcnGrSL instanceof IBMFunction_NonEggStageSTDGrowthRateSL) 
                typeGrSL = YSLStageParameters.FCN_GrSL_NonEggStageSTDGrowthRate;
            else if (fcnGrSL instanceof IBMFunction_YSL_GrowthRateSL)    
                typeGrSL = YSLStageParameters.FCN_GrSL_YSL_GrowthRate;
            
            if (fcnGrDW instanceof IBMFunction_NonEggStageSTDGrowthRateDW) 
                typeGrDW = YSLStageParameters.FCN_GrDW_NonEggStageSTDGrowthRate;
            else if (fcnGrDW instanceof IBMFunction_YSL_GrowthRateDW)    
                typeGrDW = YSLStageParameters.FCN_GrDW_YSL_GrowthRate;
            
            if (fcnVM instanceof DielVerticalMigration_FixedDepthRanges)   
                typeVM = YSLStageParameters.FCN_VM_DVM_FixedDepthRanges;
            
            if (fcnVV instanceof ConstantMovementRateFunction) 
                typeVV = YSLStageParameters.FCN_VV_YSL_Constant;
            
            if (fcnPNR instanceof IBMFunction_YSL_PNR) 
                typePNR = YSLStageParameters.FCN_PNR_YSL;
            
            if (fcnYSA instanceof IBMFunction_YSL_YSA) 
                typeYSA = YSLStageParameters.FCN_YSA_YSL;
        } else {
            //TODO: throw some error
        }
    }
    
    /*
     * Copy the values from the params map to the param variables.
     */
    private void setParameterValues() {
        isSuperIndividual = 
                params.getValue(YSLStageParameters.PARAM_isSuperIndividual,isSuperIndividual);
        horizRWP = 
                params.getValue(YSLStageParameters.PARAM_horizRWP,horizRWP);
        minStageDuration = 
                params.getValue(YSLStageParameters.PARAM_minStageDuration,minStageDuration);
        maxStageDuration = 
                params.getValue(YSLStageParameters.PARAM_maxStageDuration,maxStageDuration);
        useRandomTransitions = 
                params.getValue(YSLStageParameters.PARAM_useRandomTransitions,true);
    }
    
    /**
     *  Provides a copy of the object.  The attributes and parameters
     *  are cloned in the process, so the clone is independent of the
     *  original.
     */
    @Override
    public Object clone() {
        YSLStage clone = null;
        try {
            clone = (YSLStage) super.clone();
            clone.setAttributes(atts);//this clones atts
            clone.updateVariables();  //this sets the variables in the clone to the attribute values
            clone.setParameters(params);//this clones params
            clone.lp      = (LagrangianParticle) lp.clone();
            clone.track   = (ArrayList<Coordinate>) track.clone();
            clone.trackLL = (ArrayList<Coordinate>) trackLL.clone();
        } catch (CloneNotSupportedException ex) {
            ex.printStackTrace();
        }
        return clone;
        
    }

    /**
     *
     * @param dt - time step in seconds
     * @return
     */
    @Override
    public List<LifeStageInterface> getMetamorphosedIndividuals(double dt) {
        double dtp = 0.25*(dt/86400);//use 1/4 timestep (converted from sec to d)
        output.clear();
        List<LifeStageInterface> nLHSs;
        //SH_NEW:
        if (hasFed==true){
           if ((numTrans>0)||!isSuperIndividual){
                nLHSs = createNextLHS();
                if (nLHSs!=null) output.addAll(nLHSs);
            }
        }

        return output;
    }

    private List<LifeStageInterface> createNextLHS() {
        List<LifeStageInterface> nLHSs = null;
        try {
            //create LHS with "next" stage
            if (isSuperIndividual) {
                /** 
                 * Since this is LHS instance is a super individual, only a part 
                 * (numTrans) of it transitions to the next LHS. Thus, we need to:
                 *          1) create new LHS instance
                 *          2. assign new id to new instance
                 *          3) assign current LHS id to new LHS as parentID
                 *          4) copy current LHS origID to new LHS origID
                 *          5) set number in new LHS to numTrans for current LHS
                 *          6) reset numTrans in current LHS
                 */
                nLHSs = LHS_Factory.createNextLHSsFromSuperIndividual(typeName,this,numTrans);
                numTrans = 0.0;//reset numTrans to zero
            } else {
                /** 
                 * Since this is a single individual making a transition, we should
                 * "kill" the current LHS.  Also, the various IDs should remain
                 * the same in the new LHS since it's the same individual. Thus, 
                 * we need to:
                 *          1) create new LHS instance
                 *          2. assign current LHS id to new LHS id
                 *          3) assign current LHS parentID to new LHS parentID
                 *          4) copy current LHS origID to new LHS origID
                 *          5) kill current LHS
                 */
                nLHSs = LHS_Factory.createNextLHSsFromIndividual(typeName,this);
                alive  = false; //allow only 1 transition, so kill this stage
                active = false; //set stage inactive, also
            }
        } catch (IllegalAccessException | InstantiationException ex) {
            ex.printStackTrace();
        }
        return nLHSs;
    }
    
    @Override
    public String getReport() {
        updateAttributes();//make sure attributes are up to date
        atts.setValue(atts.PROP_track, getTrackAsString(COORDINATE_TYPE_GEOGRAPHIC));//
        return atts.getCSV();
    }

    @Override
    public String getReportHeader() {
        return atts.getCSVHeaderShortNames();
    }

    /**
     * Initializes instance variables to attribute values (via updateVariables()), 
     * then determines initial position for the lagrangian particle tracker
     * and resets the track,
     * sets horizType and vertType attributes to HORIZ_LL, VERT_H,
     * and finally calls updatePosition(), updateEnvVars(), and updateAttributes().
     */
    public void initialize() {
        updateVariables();//set instance variables to attribute values
        int hType,vType;
        hType=vType=-1;
        double xPos, yPos, zPos;
        xPos=yPos=zPos=0;
        hType      = atts.getValue(YSLStageAttributes.PROP_horizType,hType);
        vType      = atts.getValue(YSLStageAttributes.PROP_vertType,vType);
        xPos       = atts.getValue(YSLStageAttributes.PROP_horizPos1,xPos);
        yPos       = atts.getValue(YSLStageAttributes.PROP_horizPos2,yPos);
        zPos       = atts.getValue(YSLStageAttributes.PROP_vertPos,zPos);
        time       = startTime;
        numTrans   = 0.0; //set numTrans to zero
        logger.info(hType+cc+vType+cc+startTime+cc+xPos+cc+yPos+cc+zPos);
        if (i3d!=null) {
            double[] IJ = new double[] {xPos,yPos};
            if (hType==Types.HORIZ_XY) {
                IJ = i3d.getGrid().computeIJfromXY(xPos,yPos);
            } else if (hType==Types.HORIZ_LL) {
//                if (xPos<0) xPos=xPos+360;
                IJ = i3d.getGrid().computeIJfromLL(yPos,xPos);
            }
            double z = i3d.interpolateBathymetricDepth(IJ);
            logger.info("Bathymetric depth = "+z);
            double ssh = i3d.interpolateSSH(IJ);

            double K = 0;  //set K = 0 (at bottom) as default
            if (vType==Types.VERT_K) {
                if (zPos<0) {K = 0;} else
                if (zPos>i3d.getGrid().getN()) {K = i3d.getGrid().getN();} else
                K = zPos;
            } else if (vType==Types.VERT_Z) {//depths negative
                if (zPos<-z) {K = 0;} else                     //at bottom
                if (zPos>ssh) {K = i3d.getGrid().getN();} else //at surface
                K = i3d.calcKfromZ(IJ[0],IJ[1],zPos);          //at requested depth
            } else if (vType==Types.VERT_H) {//depths positive
                if (zPos>z) {K = 0;} else                       //at bottom
                if (zPos<-ssh) {K = i3d.getGrid().getN();} else //at surface
                K = i3d.calcKfromZ(IJ[0],IJ[1],-zPos);          //at requested depth
            } else if (vType==Types.VERT_DH) {//distance off bottom
                if (zPos<0) {K = 0;} else                        //at bottom
                if (zPos>z+ssh) {K = i3d.getGrid().getN();} else //at surface
                K = i3d.calcKfromZ(IJ[0],IJ[1],-(z-zPos));       //at requested distance off bottom
            }
            lp.setIJK(IJ[0],IJ[1],K);
            //reset track array
            track.clear();
            trackLL.clear();
            //set horizType to lat/lon and vertType to depth
            atts.setValue(YSLStageAttributes.PROP_horizType,Types.HORIZ_LL);
            atts.setValue(YSLStageAttributes.PROP_vertType,Types.VERT_H);
            //interpolate initial position and environmental variables
            double[] pos = lp.getIJK();
            updatePosition(pos);
            interpolateEnvVars(pos);
            updateAttributes(); 
        }
    }
    
    @Override
    public void step(double dt) throws ArrayIndexOutOfBoundsException {
        //WTS_NEW 2012-07-26:{
        double[] pos = lp.getIJK();
        double T0 = i3d.interpolateTemperature(pos);
        
      //SH-Prey Stuff  
        copepods    = i3d.interpolateValue(pos,Cop,Interpolator3D.INTERP_VAL);
        euphausiids = i3d.interpolateValue(pos,Eup,Interpolator3D.INTERP_VAL);
        neocalanus  = i3d.interpolateValue(pos,NCa,Interpolator3D.INTERP_VAL);
               
        double[] res = calcW(pos,dt);//calc w and attached indicator
        double w     = Math.signum(dt)*res[0];
        attached     = res[1]<0;
        if (attached) pos[2] = 0;//set individual on bottom
        double[] uv  = calcUV(pos,dt);//calculate orizontal movement components
        if (attached){
            lp.setIJK(pos[0], pos[1], pos[2]);
        } else {
            //}:WTS_NEW 2012-07-26
            //do lagrangian particle tracking
            lp.setU(uv[0],lp.getN());
            lp.setV(uv[1],lp.getN());
            lp.setW(w,    lp.getN());
            //now do predictor step
            lp.doPredictorStep();
            //assume same daytime status, but recalc depth and revise W 
            pos = lp.getPredictedIJK();
            depth = -i3d.calcZfromK(pos[0],pos[1],pos[2]);
            if (debug) logger.info("Depth after predictor step = "+depth);
            //w = calcW(dt,lp.getNP1())+r; //set swimming rate for predicted position
            lp.setU(uv[0],lp.getNP1());
            lp.setV(uv[1],lp.getNP1());
            lp.setW(w,    lp.getNP1());
            //now do corrector step
            lp.doCorrectorStep();
            pos = lp.getIJK();
            if (debug) logger.info("Depth after corrector step = "+(-i3d.calcZfromK(pos[0],pos[1],pos[2])));
        }
        
        time += dt;
        double dtday = dt/86400;//bio model timestep in days
        
        //get effective temperature as average temp at new and old locations
        double T1 = i3d.interpolateTemperature(pos);
        double T = 0.5 * (T0 + T1);
        if(T<=0.0) T=0.01; 

        //Days to 100% mortality        
        durPNR  = (Double) fcnPNR.calculate(T);//only 1 alternative function currently defined
        progPNR += dtday/durPNR;//integrated criterion for point-of-no return (progPNR=1)
        
        if (progPNR>=1.0){
            alive=false;//larva passes PNR without feeding and dies of starvation
            active=false;
        } else {
            //Days to YSA (when it is ready to feed)
            durYSA  = (Double) fcnYSA.calculate(T);//only 1 alternative function currently defined        
            if (progYSA<1.0) progYSA += dtday/durYSA;//integrated criterion for yolk-sac absorption (progYSA=1)
            if ((ageYSA<0)&&(progYSA>=1.0))
                ageYSA = ageInStage; //Age at which feeding is possible

            //growth is same for feeding via ysa or active feeding 
            if (typeGrSL==YSLStageParameters.FCN_GrSL_YSL_GrowthRate)
                grSL = (Double) fcnGrSL.calculate(T);
            else if (typeGrSL==YSLStageParameters.FCN_GrSL_NonEggStageSTDGrowthRate)
                grSL = (Double) fcnGrSL.calculate(new Double[]{T,std_len});

            if (typeGrDW==YSLStageParameters.FCN_GrDW_YSL_GrowthRate)
                grDW = (Double) fcnGrDW.calculate(T);
            else if (typeGrDW==YSLStageParameters.FCN_GrDW_NonEggStageSTDGrowthRate)
                grDW = (Double) fcnGrDW.calculate(new Double[]{T,dry_wgt});

            if (progYSA<1.0){
                //yolk-sac absorption is incomplete
                std_len += grSL*dtday;            //mm
                dry_wgt *= Math.exp(grDW * dtday);//mg
            } else {
                //yolk sac absorption is complete
                //calculate whether or not first feeding occurs                
                if (useFirstFeedingSH){
                    //SH approach
                    //This approach assumes prFeed (the Lifetime Distribution Function) 
                    //increases linearly from 0 to 1 as time/age increases from YSA to PNR,
                    //but this occurs necessarily only in the case of constant temperature
                    //(in which case durPNR-durYSA is a constant).
                    prFeed  += dtday/(durPNR-durYSA);
                    prNotFed = 1.0-prFeed;
                    double b = Math.random();
                    logger.info("Check on first feeding for id "+id+": "+rndFeed+" <= "+prFeed+"?");
                    if (rndFeed<=prFeed) hasFed = true;//feeding occurs, will transition to FDL stage
                    //growth occurs regardless of feeding (seems unrealistic)
                    std_len += grSL*dtday;
                    dry_wgt *= Math.exp(grDW * dtday);
                } else {
                    //WTS approach based on survival/failure analysis
                    //hazard function for first feeding is based on prey encounter rate (search volume x abundance density)
                    double svr  = Math.PI*(std_len*std_len*(1.0e-6))*Math.abs(w);//search volume rate (m^3/s)
                    double fHF  = svr*(copepods/indivCopWgt);//instantaneous feeding hazard rate (1/s)
                    fCumHazFcn += fHF*dt;                   //cumulative hazard function for first feeding (note: dt, not dtday)
                    prNotFed = Math.exp(-fCumHazFcn);
                    if (rndFeed>prNotFed) {
                        //feeding occurs
                        hasFed = true;//will transition to FDL stage
                        std_len += grSL*dtday;
                        dry_wgt *= Math.exp(grDW * dtday);
                    } else {
                        //growth does not occur if ysa is completed but feeding has not begun
                    }
                }
            }
        }
          
        updateNum(dt);
        updateAge(dt);
        updatePosition(pos);
        interpolateEnvVars(pos);
        //check for exiting grid
        if (i3d.isAtGridEdge(pos,tolGridEdge)){
            alive=false;
            active=false;
        }
        if (debug) {
            logger.info(toString());
        }
        updateAttributes(); //update the attributes object w/ nmodified values
    }
    
    /**
     * Function to calculate vertical movement rate (m/s).
     * 
     * @param pos - position vector
     * @param dt - time step
     * 
     * @return double[] with
     *  [0] - w, individual active vertical movement velocity (m/s)
     *  [1] - flag indicating whether individual is attached to bottom(< 0) or not (>0)
     */
    private double[] calcW(double[] pos, double dt){
        //compute vertical velocity
        double[] res = null;
        double w = 0;
        if (typeVM==YSLStageParameters.FCN_VM_DVM_FixedDepthRanges) {
            //fcnVM instanceof wts.models.DisMELS.IBMFunctions.Movement.DielVerticalMigration_FixedDepthRanges
            //calculate the vertical movement rate
            if (typeVV==YSLStageParameters.FCN_VV_YSL_Constant) {
                //fcnVV instanceof wts.models.DisMELS.IBMFunctions.SwimmingBehavior.ConstantMovementRateFunction
                /**
                * @param vars - double[]{dt}.
                * @return     - movement rate as a Double 
                */
                w = (Double) fcnVV.calculate(new double[]{dt});
            }
            /**
            * Compute time of local sunrise, sunset and solar noon (in minutes, UTC) 
            * for given lon, lat, and time (in Julian day-of-year).
            *@param lon : longitude of position (deg Greenwich, prime meridian)
            *@param lat : latitude of position (deg)
            *@param time : day-of-year (1-366, fractional part indicates time-of-day)
            *@return double[5] = [0] time of sunrise (min UTC from midnight)
            *                    [1] time of sunset (min UTC from midnight)
            *                    [2] time of solarnoon (min UTC from midnight)
            *                    [3] solar declination angle (deg)
            *                    [4] solar zenith angle (deg)
            * If sunrise/sunset=NaN then its either 24-hr day or night 
            * (if lat*declination>0, it's summer in the hemisphere, hence daytime). 
            * Alternatively, if the solar zenith angle > 90.833 deg, then it is night.
            */
            CalendarIF cal = null;
            double[] ss = null;
            try {
                cal = GlobalInfo.getInstance().getCalendar();
                ss = DateTimeFunctions.computeSunriseSunset(lon,lat,cal.getYearDay());
            } catch(java.lang.NullPointerException ex){
                logger.info("NullPointerException for EggStage id: "+id);
                logger.info("lon: "+lon+". lat: "+lat+". yearday: "+cal.getYearDay());
                logger.info(ex.getMessage());
            }
            /**
            * @param vars - the inputs variables as a double[] array with elements
            *                  dt          - [0] - integration time step
            *                  depth       - [1] - current depth of individual
            *                  total depth - [2] - total depth at location
            *                  w           - [3] - active vertical swimming speed outside preferred depth range
            *                  lightLevel  - [4] - value >= 0 indicates daytime, otherwise night 
            * @return     - double[] with elements
            *              w        - individual active vertical movement velocity
            *              attached - flag indicating whether individual is attached to bottom(< 0) or not (>0)
            */
            double td = i3d.interpolateBathymetricDepth(lp.getIJK());
            res = (double[]) fcnVM.calculate(new double[]{dt,depth,td,w,90.833-ss[4]});
        }
        return res;
    }
    
    /**
     * Function to calculate horizontal movement rates.
     * 
     * @param pos - position vector
     * @param dt - time step
     * @return 
     */
    public double[] calcUV(double[] pos, double dt) {
        //calculate horizontal movement
        double[] uv = {0.0,0.0};
        if (!attached){
            if ((horizRWP>0)&&(Math.abs(dt)>0)) {
                double r = Math.sqrt(horizRWP/Math.abs(dt));
                uv[0] += r*rng.computeNormalVariate(); //stochastic swimming rate
                uv[1] += r*rng.computeNormalVariate(); //stochastic swimming rate
                if (debug) System.out.print("uv: "+r+"; "+uv[0]+", "+uv[1]+"\n");
            }
        }
        //return the result
        return new double[]{Math.signum(dt)*uv[0],Math.signum(dt)*uv[1]};
    }
    //WTS_NEW 2012-07-26:{

    /**
     *
     * @param dt - time step in seconds
     */
    private void updateAge(double dt) {
        age        += dt/86400;
        ageInStage += dt/86400;
        if (ageInStage>maxStageDuration) {
            alive = false;
            active = false;
        }
    }

    /**
     *
     * @param dt - time step in seconds
     */
    private void updateNum(double dt) {
        //{WTS_NEW 2012-07-26:
        double mortalityRate = 0.0D;//in unis of [days]^-1
        if (typeMort==YSLStageParameters.FCN_Mortality_ConstantMortalityRate){
            //fcnMortality instanceof ConstantMortalityRate
            mortalityRate = (Double)fcnMortality.calculate(null);
        } else 
        if (typeMort==YSLStageParameters.FCN_Mortality_InversePowerLawMortalityRate){
            //fcnMortality instanceof InversePowerLawMortalityRate
            mortalityRate = (Double)fcnMortality.calculate(std_len);//using std_len as covariate for mortality
        }
        double totRate = mortalityRate;
        if ((ageInStage>=minStageDuration)) {
            totRate += stageTransRate;
            //apply mortality rate to previous number transitioning and
            //add in new transitioners
            numTrans = numTrans*Math.exp(-dt*mortalityRate/86400)+
                    (stageTransRate/totRate)*number*(1-Math.exp(-dt*totRate/86400));
        }
        number = number*Math.exp(-dt*totRate/86400);
        //}: WTS_NEW 2012-07-26
    }
    
    private void updatePosition(double[] pos) {
        bathym     =  i3d.interpolateBathymetricDepth(pos);
        depth      = -i3d.calcZfromK(pos[0],pos[1],pos[2]);
        lat        =  i3d.interpolateLat(pos);
        lon        =  i3d.interpolateLon(pos);
        gridCellID = ""+Math.round(pos[0])+"_"+Math.round(pos[1]);
        updateTrack();
    }
    
    private void interpolateEnvVars(double[] pos) {
        temperature = i3d.interpolateTemperature(pos);
        salinity    = i3d.interpolateSalinity(pos);
        if (i3d.getPhysicalEnvironment().getField("rho")!=null) 
            rho  = i3d.interpolateValue(pos,"rho");
        else rho = 0.0;
    }

    @Override
    public double getStartTime() {
        return startTime;
    }

    @Override
    public void setStartTime(double newTime) {
        startTime = newTime;
        time      = startTime;
        atts.setValue(YSLStageAttributes.PROP_startTime,startTime);
        atts.setValue(YSLStageAttributes.PROP_time,time);
    }

    @Override
    public boolean isActive() {
        return active;
    }

    @Override
    public void setActive(boolean b) {
        active = b;
        atts.setActive(b);
    }

    @Override
    public boolean isAlive() {
        return alive;
    }

    @Override
    public void setAlive(boolean b) {
        alive = b;
        atts.setAlive(b);
    }

    @Override
    public String getAttributesClassName() {
        return attributesClass;
    }

    @Override
    public String getParametersClassName() {
        return parametersClass;
    }

    @Override
    public String[] getNextLHSClassNames() {
        return nextLHSClasses;
    }

    @Override
    public String getPointFeatureTypeClassName() {
        return pointFTClass;
    }

    @Override
    public String[] getSpawnedLHSClassNames() {
        return spawnedLHSClasses;
    }

    @Override
    public List<LifeStageInterface> getSpawnedIndividuals() {
        output.clear();
        return output;
    }

    @Override
    public boolean isSuperIndividual() {
        return isSuperIndividual;
    }
    
    /**
     * Updates attribute values defined for this abstract class. 
     */
    @Override
    protected void updateAttributes() {
        super.updateAttributes();
        atts.setValue(YSLStageAttributes.PROP_attached,attached);
        atts.setValue(YSLStageAttributes.PROP_SL,std_len);
        atts.setValue(YSLStageAttributes.PROP_DW,dry_wgt);
        atts.setValue(YSLStageAttributes.PROP_grSL,grSL);
        atts.setValue(YSLStageAttributes.PROP_grDW,grDW);
        atts.setValue(YSLStageAttributes.PROP_temperature,temperature);
        atts.setValue(YSLStageAttributes.PROP_salinity,salinity);
        atts.setValue(YSLStageAttributes.PROP_rho,rho);
        atts.setValue(YSLStageAttributes.PROP_copepod,copepods);
        atts.setValue(YSLStageAttributes.PROP_euphausiid,euphausiids);
        atts.setValue(YSLStageAttributes.PROP_neocalanus,neocalanus);
        atts.setValue(YSLStageAttributes.PROP_progYSA,progYSA);
        atts.setValue(YSLStageAttributes.PROP_progPNR,progPNR);
        atts.setValue(YSLStageAttributes.PROP_prNotFed,prNotFed);
    }

    /**
     * Updates local variables from the attributes.  
     */
    @Override
    protected void updateVariables() {
        super.updateVariables();
        attached    = atts.getValue(YSLStageAttributes.PROP_attached,attached);
        std_len     = atts.getValue(YSLStageAttributes.PROP_SL,std_len); 
        dry_wgt     = atts.getValue(YSLStageAttributes.PROP_DW,dry_wgt); 
        grSL        = atts.getValue(YSLStageAttributes.PROP_grSL,grSL); 
        grDW        = atts.getValue(YSLStageAttributes.PROP_grDW,grDW); 
        temperature = atts.getValue(YSLStageAttributes.PROP_temperature,temperature);
        salinity    = atts.getValue(YSLStageAttributes.PROP_salinity,salinity);
        rho         = atts.getValue(YSLStageAttributes.PROP_rho,rho);
        copepods    = atts.getValue(YSLStageAttributes.PROP_copepod,copepods);
        euphausiids = atts.getValue(YSLStageAttributes.PROP_euphausiid,euphausiids);
        neocalanus  = atts.getValue(YSLStageAttributes.PROP_neocalanus,neocalanus);
        progYSA     = atts.getValue(YSLStageAttributes.PROP_progYSA,progYSA); 
        progPNR     = atts.getValue(YSLStageAttributes.PROP_progPNR,progPNR); 
        prNotFed    = atts.getValue(YSLStageAttributes.PROP_prNotFed,prNotFed); 
    }

}
