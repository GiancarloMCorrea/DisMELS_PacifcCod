/**
 * FDLStage.java
 *
 * Revisions:
 * 20181011:  1. Corrected gL formula
 *            2. Removed fcnDevelopment to match Parameters class.
 *            3. Removed fcnVV because calculation for w is hard-wired in calcUVW.
 *            4. Added "attached" as new attribute (necessary with updated DisMELS).
 *            5. Removed "diam" since it's replaced by "length"
 * 20190722: 1. Removed fields associated with egg stage attributes "devStage" and "density"
 *
 */

package sh.pcod.FDLStage;

import com.vividsolutions.jts.geom.Coordinate;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.List;
import org.openide.util.lookup.ServiceProvider;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageBIOENGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateSL;
import wts.models.DisMELS.IBMFunctions.Mortality.ConstantMortalityRate;
import wts.models.DisMELS.IBMFunctions.Mortality.InversePowerLawMortalityRate;
import wts.models.DisMELS.framework.*;
import wts.models.DisMELS.framework.IBMFunctions.IBMFunctionInterface;
import wts.models.utilities.DateTimeFunctions;
import wts.roms.model.LagrangianParticle;
import sh.pcod.YSLStage.YSLStageAttributes;
import wts.models.DisMELS.IBMFunctions.Movement.DielVerticalMigration_FixedDepthRanges;
import wts.models.utilities.CalendarIF;
import wts.roms.model.Interpolator3D;
// import wts.roms.model.Interpolator2D;

/**
 *
 * @author William Stockhausen
 */
@ServiceProvider(service=LifeStageInterface.class)
public class FDLStage extends AbstractLHS {
    
        //Static fields    
            //  Static fields new to this class
    /* flag to do debug operations */
    public static boolean debug = false;
    /* Class for attributes SH_NEW */
    public static final String attributesClass = 
            sh.pcod.FDLStage.FDLStageAttributes.class.getName();
    /* Class for parameters */
    public static final String parametersClass = 
            sh.pcod.FDLStage.FDLStageParameters.class.getName();
    /* Class for feature type for point positions */
    public static final String pointFTClass = 
            wts.models.DisMELS.framework.LHSPointFeatureType.class.getName();
    /* Classes for next LHS SH_NEW Add fdl*/
    public static final String[] nextLHSClasses = new String[]{ 
            sh.pcod.FDLpfStage.FDLpfStage.class.getName()};
    /* Classes for spawned LHS */
    public static final String[] spawnedLHSClasses = new String[]{};
    
    /* string identifying environmental field with copepod densities */
    private static final String Cop = "Cop";
    /* string identifying environmental field with euphausiid densities */
    private static final String Eup = "Eup";
    /* string identifying environmental field with neocalanus densities */
    private static final String NCa = "NCa";
    /* string identifying environmental field with large phytoplankton densities */
    private static final String PhL = "PhL";
    /* string identifying environmental field with small phytoplankton densities */
    private static final String PhS = "PhS";
    /* string identifying environmental field with u surface stress */
    private static final String Su = "Su";
    /* string identifying environmental field with v surface stress */
    private static final String Sv = "Sv";

    //Instance fields
            //  Fields hiding ones from superclass
    /* life stage atrbutes object */
    protected FDLStageAttributes atts = null;
    /* life stage parameters object */
    protected FDLStageParameters params = null;
    
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
    /** stomach state (units) */
    protected double stmsta = 0;
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
    /** in situ small copepod density (mg/m^3, dry wt) */
     protected double copepod; 
     /** in situ euphausiid density (mg/m^3, dry wt) */
    protected double euphausiid = 0;
     /** in situ neocalanoid density (mg/m^3, dry wt) */
    protected double neocalanus = 0;
     /** in situ large phytoplankton density (mg/m^3) */
    // protected double phytoL = 0;
     /** in situ small phytoplankton density (mg/m^3) */
    // protected double phytoS = 0;
     /** in situ u surface stress (N/m^2) */
    // protected double tauX = 0;
     /** in situ v surface stress (N/m^2) */
    // protected double tauY = 0;

            //other fields
    /** number of individuals transitioning to next stage */
    private double numTrans;  
    /** FDL Size at flexion */
    protected double flexion=13.5;
    /** in situ temperature */
    double T;
    
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
    
    private int typeMort = 0;//integer indicating mortality function
    private int typeGrSL = 0;//integer indicating SL growth function
    private int typeGrDW = 0;//integer indicating DW growth function
    private int typeVM   = 0;//integer indicating vertical movement function
    private int typeVV   = 0;//integer indicating vertical velocity function

    private static final Logger logger = Logger.getLogger(FDLStage.class.getName());
    
    /**
     * Creates a new instance of FDLStage.  
     *  This constructor should be used ONLY to obtain
     *  the class names of the associated classes.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     */
    public FDLStage() {
        super("");
        super.atts = atts;
        super.params = params;
    }
    
    /**
     * Creates a new instance of FDLStage with the given typeName.
     * A new id number is calculated in the superclass and assigned to
     * the new instance's id, parentID, and origID. 
     * 
     * The attributes are vanilla.  Initial attribute values should be set,
     * then initialize() should be called to initialize all instance variables.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     * @param typeName
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     */
    public FDLStage(String typeName) 
                throws InstantiationException, IllegalAccessException {
        super(typeName);
        atts = new FDLStageAttributes(typeName);
        atts.setValue(FDLStageAttributes.PROP_id,id);
        atts.setValue(FDLStageAttributes.PROP_parentID,id);
        atts.setValue(FDLStageAttributes.PROP_origID,id);
        params = (FDLStageParameters) LHS_Factory.createParameters(typeName);
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
    public FDLStage createInstance(String[] strv) 
                        throws InstantiationException, IllegalAccessException {
        LifeStageAttributesInterface theAtts = LHS_Factory.createAttributes(strv);
        FDLStage lhs = createInstance(theAtts);
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
    public FDLStage createInstance(LifeStageAttributesInterface theAtts)
                        throws InstantiationException, IllegalAccessException {
        FDLStage lhs = null;
        if (theAtts instanceof FDLStageAttributes) {
            lhs = new FDLStage(theAtts.getTypeName());
            long newID = lhs.id;//save id of new instance
            lhs.setAttributes(theAtts);
            if (lhs.atts.getID()==-1) {
                //constructing new individual, so reset id values to those of new
                lhs.id = newID;
                lhs.atts.setValue(FDLStageAttributes.PROP_id,newID);
            }
            newID = (Long) lhs.atts.getValue(FDLStageAttributes.PROP_parentID);
            if (newID==-1) {
                lhs.atts.setValue(FDLStageAttributes.PROP_parentID,newID);
            }
            newID = (Long) lhs.atts.getValue(FDLStageAttributes.PROP_origID);
            if (newID==-1) {
                lhs.atts.setValue(FDLStageAttributes.PROP_origID,newID);
            }
        }
        lhs.initialize();//initialize instance variables
        return lhs;
    }

    /**
     *  Returns the associated attributes.  
     */
    @Override
    public FDLStageAttributes getAttributes() {
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
        aid = atts.getValue(FDLStageAttributes.PROP_id, id);
        if (aid==-1) {
            //change atts id to lhs id
            atts.setValue(FDLStageAttributes.PROP_id, id);
        } else {
            //change lhs id to atts id
            id = aid;
        }
        aid = atts.getValue(FDLStageAttributes.PROP_parentID, id);
        if (aid==-1) {
            atts.setValue(FDLStageAttributes.PROP_parentID, id);
        }
        aid = atts.getValue(FDLStageAttributes.PROP_origID, id);
        if (aid==-1) {
            atts.setValue(FDLStageAttributes.PROP_origID, id);
        }
        initialize();//initialize instance variables
    }

    /**
     * Sets the attributes for the instance by copying values from the input.
     * This does NOT change the typeName of the LHS instance (or the associated 
     * LHSAttributes instance) on which the method is called.
     * Note that ALL attributes are copied, so id, parentID, and origID are copied
     * as well. 
     *  Side effects:
     *      updateVariables() is called to update instance variables.
     *      Instance field "id" is also updated.
     * @param newAtts - should be instance of SimplePelagicLHSAttributes
     */
    @Override
    public void setAttributes(LifeStageAttributesInterface newAtts) {
        if (newAtts instanceof FDLStageAttributes) {
            FDLStageAttributes oldAtts = (FDLStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));}
        else if(newAtts instanceof YSLStageAttributes) {
            YSLStageAttributes oldAtts = (YSLStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));
            //all FDL attributes are also YSL attributes, so no need to do anything further
        } else {
            //TODO: should throw an error here
            logger.info("setAttributes(): no match for attributes type:"+newAtts.toString());
        }
        id = atts.getValue(FDLStageAttributes.PROP_id, id);
        //updateVariables();//this gets done in setInfoFrom... after call to this method, so no need to do it here
    }
    
    /**
     *  Sets the associated attributes object. Use this after creating an LHS instance
     * as an "output" from another LHS that is functioning as an ordinary individual.
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
        atts.setValue(FDLStageAttributes.PROP_ageInStage, 0.0);//reset age in stage
        atts.setValue(FDLStageAttributes.PROP_active,true);    //set active to true
        atts.setValue(FDLStageAttributes.PROP_alive,true);     //set alive to true
        id = atts.getID(); //reset id for current LHS to one from old LHS

        //update local variables to capture changes made here
        updateVariables();
    }
    
    /**
     *  Sets the associated attributes object. Use this after creating an LHS instance
     * as an "output" from another LHS that is functioning as a super individual.
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
        atts.setValue(FDLStageAttributes.PROP_id,idc);//reset id to one for current LHS
        atts.setValue(FDLStageAttributes.PROP_parentID,
                      oldAtts.getValue(LifeStageAttributesInterface.PROP_id));//copy old id to parentID
        atts.setValue(FDLStageAttributes.PROP_number, numTrans);//set number to numTrans
        atts.setValue(FDLStageAttributes.PROP_ageInStage, 0.0); //reset age in stage
        atts.setValue(FDLStageAttributes.PROP_active,true);     //set active to true
        atts.setValue(FDLStageAttributes.PROP_alive,true);      //set alive to true
            
        //update local variables to capture changes made here
        updateVariables();
    }

    /**
     *  Returns the associated parameters.  
     * @return 
     */
    @Override
    public FDLStageParameters getParameters() {
        return params;
    }

    /**
     * Sets the parameters for the instance to a cloned version of the input.
     * @param newParams - should be instance of FDLStageParameters
     */
    @Override
    public void setParameters(LifeStageParametersInterface newParams) {
        if (newParams instanceof FDLStageParameters) {
            params = (FDLStageParameters) newParams;
            super.params = params;
            setParameterValues();
            fcnMortality = params.getSelectedIBMFunctionForCategory(FDLStageParameters.FCAT_Mortality);
            fcnGrSL = params.getSelectedIBMFunctionForCategory(FDLStageParameters.FCAT_GrowthSL);
            fcnGrDW = params.getSelectedIBMFunctionForCategory(FDLStageParameters.FCAT_GrowthDW);
            fcnVM   = params.getSelectedIBMFunctionForCategory(FDLStageParameters.FCAT_VerticalMovement);
            fcnVV   = params.getSelectedIBMFunctionForCategory(FDLStageParameters.FCAT_VerticalVelocity);
            
            if (fcnMortality instanceof ConstantMortalityRate)
                typeMort = FDLStageParameters.FCN_Mortality_ConstantMortalityRate;
            else if (fcnMortality instanceof InversePowerLawMortalityRate)
                typeMort = FDLStageParameters.FCN_Mortality_InversePowerLawMortalityRate;
            
            if (fcnGrSL instanceof IBMFunction_NonEggStageSTDGrowthRateSL) 
                typeGrSL = FDLStageParameters.FCN_GrSL_NonEggStageSTDGrowthRate;
            else if (fcnGrSL instanceof IBMFunction_FDL_GrowthRateSL)    
                typeGrSL = FDLStageParameters.FCN_GrSL_FDL_GrowthRate;
            
            if (fcnGrDW instanceof IBMFunction_NonEggStageSTDGrowthRateDW) {
                typeGrDW = FDLStageParameters.FCN_GrDW_NonEggStageSTDGrowthRate;
            } else if (fcnGrDW instanceof IBMFunction_FDL_GrowthRateDW) {  
                typeGrDW = FDLStageParameters.FCN_GrDW_FDL_GrowthRate;
            } else if (fcnGrDW instanceof IBMFunction_NonEggStageBIOENGrowthRateDW) {
                typeGrDW = FDLStageParameters.FCN_GrDW_NonEggStageBIOENGrowthRate;
            }
            
            if (fcnVM instanceof DielVerticalMigration_FixedDepthRanges)   
                typeVM = FDLStageParameters.FCN_VM_DVM_FixedDepthRanges;
            
            if (fcnVV instanceof IBMFunction_FDL_VerticalSwimmingSpeed) 
                typeVV = FDLStageParameters.FCN_VV_FDL_VerticalSwimmingSpeed;
        } else {
            //TODO: throw some error
        }
    }
    
    /*
     * Copy the values from the params map to the param variables.
     */
    private void setParameterValues() {
        isSuperIndividual = 
                params.getValue(FDLStageParameters.PARAM_isSuperIndividual,isSuperIndividual);
        horizRWP = 
                params.getValue(FDLStageParameters.PARAM_horizRWP,horizRWP);
        minStageDuration = 
                params.getValue(FDLStageParameters.PARAM_minStageDuration,minStageDuration);
        maxStageDuration = 
                params.getValue(FDLStageParameters.PARAM_maxStageDuration,maxStageDuration);
        useRandomTransitions = 
                params.getValue(FDLStageParameters.PARAM_useRandomTransitions,true);
    }
    
    /**
     *  Provides a copy of the object.  The attributes and parameters
     *  are cloned in the process, so the clone is independent of the
     *  original.
     */
    @Override
    public Object clone() {
        FDLStage clone = null;
        try {
            clone = (FDLStage) super.clone();
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
        
        if(std_len >= flexion) {
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
//        atts.setValue(SimplePelagicLHSAttributes.PARAM_id,id);//TODO: should do this beforehand!!
        updateVariables();//set instance variables to attribute values
        int hType,vType;
        hType=vType=-1;
        double xPos, yPos, zPos;
        xPos=yPos=zPos=0;
        hType      = atts.getValue(FDLStageAttributes.PROP_horizType,hType);
        vType      = atts.getValue(FDLStageAttributes.PROP_vertType,vType);
        xPos       = atts.getValue(FDLStageAttributes.PROP_horizPos1,xPos);
        yPos       = atts.getValue(FDLStageAttributes.PROP_horizPos2,yPos);
        zPos       = atts.getValue(FDLStageAttributes.PROP_vertPos,zPos);
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
            atts.setValue(FDLStageAttributes.PROP_horizType,Types.HORIZ_LL);
            atts.setValue(FDLStageAttributes.PROP_vertType,Types.VERT_H);
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
        double[] pos2d = new double[2];
        pos2d[0] = pos[0];
        pos2d[1] = pos[1];
        //System.out.print("uv: "+r+"; "+uv[0]+", "+uv[1]+"\n");
        T = i3d.interpolateTemperature(pos);
        if(T<=0.0) T=0.01; 
             //SH-Prey Stuff  
        copepod    = i3d.interpolateValue(pos,Cop,Interpolator3D.INTERP_VAL);
        //euphausiid = i3d.interpolateValue(pos,Eup,Interpolator3D.INTERP_VAL);
        neocalanus = i3d.interpolateValue(pos,NCa,Interpolator3D.INTERP_VAL);
        // ADD HERE OTHER ZOOPLANKTON PREY ITEMS
        double phytoL = i3d.interpolateValue(pos,PhL,Interpolator3D.INTERP_VAL);
        double phytoS = i3d.interpolateValue(pos,PhS,Interpolator3D.INTERP_VAL);
        double tauX = i3d.interpolateValue(pos2d,Su,Interpolator3D.INTERP_VAL); // 3D interpolator but should use 2D internally
        double tauY = i3d.interpolateValue(pos2d,Sv,Interpolator3D.INTERP_VAL); // 3D interpolator but should use 2D internally
        // double tauX = -0.15; // TODO: link with ROMS output. Surface stress In N/m^2
        // double tauY = 0.1; // TODO: link with ROMS output. Surface stress In N/m^2
        double chlorophyll = (phytoL/25) + (phytoS/65); // calculate chlorophyll (mg/m^-3) 
        // values 25 and 65 based on Kearney et al 2018 Table A4
        //double bathy = i3d.interpolateBathymetricDepth(pos);

        double[] uvw = calcUVW(pos,dt);//this also sets "attached" and may change pos[2] to 0
        euphausiid = uvw[2]; // print out this value
        //PRINT UVW
        if (attached){
            lp.setIJK(pos[0], pos[1], pos[2]);
        } else {
            //}:WTS_NEW 2012-07-26
            //do lagrangian particle tracking
            lp.setU(uvw[0],lp.getN());
            lp.setV(uvw[1],lp.getN());
            lp.setW(uvw[2],lp.getN());
            //now do predictor step
            lp.doPredictorStep();
            //assume same daytime status, but recalc depth and revise W 
            pos = lp.getPredictedIJK();
            //PRINT HERE
            depth = -i3d.calcZfromK(pos[0],pos[1],pos[2]);
            if (debug) logger.info("Depth after predictor step = "+depth);
            //w = calcW(dt,lp.getNP1())+r; //set swimming rate for predicted position
            lp.setU(uvw[0],lp.getNP1());
            lp.setV(uvw[1],lp.getNP1());
            lp.setW(uvw[2],lp.getNP1());
            //now do corrector step
            lp.doCorrectorStep();
            pos = lp.getIJK();
            if (debug) logger.info("Depth after corrector step = "+(-i3d.calcZfromK(pos[0],pos[1],pos[2])));
        }
        
        // BEGIN OF BIOEN CHANGES --------------------------------------------------------------

        time += dt;
        double dtday = dt/86400; //time setp in days

        // Create objects for output of BIOEN:
        Double[] bioEN_output = null; // output object of BIOEN calculation
        double gr_mg_fac = 0; // factor to avoid dry_wgt in IF 
        double meta = 0;
        double sum_ing = 0;
        double assi = 0;
        double old_dry_wgt = dry_wgt; // save previous dry_wgt
        // double lat = pos[2]; // lat value not in this way
        double old_std_len = std_len;

        // Length:
        if (typeGrSL==FDLStageParameters.FCN_GrSL_NonEggStageSTDGrowthRate) {
            grSL = (Double) fcnGrSL.calculate(new Double[]{T,std_len});
            std_len += grSL*dtday;
        }
        else if (typeGrSL==FDLStageParameters.FCN_GrSL_FDL_GrowthRate) {
            grSL = (Double) fcnGrSL.calculate(T);
            std_len += grSL*dtday;
        }

        // Weight:
        if (typeGrDW==FDLStageParameters.FCN_GrDW_NonEggStageSTDGrowthRate) {
            grDW = (Double) fcnGrDW.calculate(new Double[]{T,dry_wgt});
            gr_mg_fac = dry_wgt*(Math.exp(grDW*dtday) - 1);
            dry_wgt += gr_mg_fac;
        }
        if (typeGrDW==FDLStageParameters.FCN_GrDW_FDL_GrowthRate) {
            grDW = (Double) fcnGrDW.calculate(T);
            gr_mg_fac = dry_wgt*(Math.exp(grDW*dtday) - 1);
            dry_wgt += gr_mg_fac;
        }
        if (typeGrDW==FDLStageParameters.FCN_GrDW_NonEggStageBIOENGrowthRate){

            // Light (begin):
            // create object for light calculation:
            double[] eb2 = new double[2]; // K parameter and second part of Eb equation
            double eb = 0; // create Eb object
            CalendarIF cal = null;
            double[] ltemp = new double[3];
            double[] ltemp2 = new double[2];
            cal = GlobalInfo.getInstance().getCalendar(); // to calculate julian day
            ltemp = IBMFunction_NonEggStageBIOENGrowthRateDW.calcLightQSW(lat,cal.getYearDay()); // see line 713 in ibm.py
            double maxLight = ltemp[0]/0.217; // see line 714 in ibm.py
            ltemp2 = IBMFunction_NonEggStageBIOENGrowthRateDW.calcLightSurlig(lat,cal.getYearDay(), maxLight); // see line 715 in ibm.py
            eb2 = IBMFunction_NonEggStageBIOENGrowthRateDW.calcLight(chlorophyll, depth, bathym); // second part of Eb equation
            eb = 0.42*ltemp2[1]*eb2[1]; // see line 727 in ibm.py. This is Eb. 0.42 as in Kearney et al 2020 Eq A14
            // Light (end):

            // Turbulence and wind (begin)
            double windX = Math.abs(Math.sqrt(Math.abs(tauX)/(1.3*1.2E-3))); // Wind velocity In m/s
            double windY = Math.abs(Math.sqrt(Math.abs(tauY)/(1.3*1.2E-3))); // Wind velocity In m/s
            // Turbulence and wind (end)

            // Bioenergetic growth calculation:
            bioEN_output = (Double[]) fcnGrDW.calculate(new Double[]{T,dry_wgt,dt,dtday,old_std_len,eb,windX,windY,depth,stmsta,copepod,eb2[0]}); //should length be at t or t-1?
            grDW = bioEN_output[0]; // grDW is gr_mg in TROND here
            meta = bioEN_output[1];
            sum_ing = bioEN_output[2];
            assi = bioEN_output[3];
            double costRateOfMetabolism = 0.5; // check this number
            double activityCost = 0.5*meta*costRateOfMetabolism; // TODO: (diffZ/maxDiffZ) = 0.5, but this should change based on vertical movement

            // Update values:
            stmsta = Math.max(0, Math.min(0.06*dry_wgt, stmsta + sum_ing)); // gut_size= 0.06. TODO: check if sum_ing is 500 approx makes sense
            gr_mg_fac = Math.min(grDW + meta, stmsta*assi) - meta - activityCost; // Here grDW is as gr_mg in TROND
            
            // TODO: decide how many days a larva can survive without food

            // TODO: try paper formulation of what growth equations to use.

            // new weight in mg:
            dry_wgt += gr_mg_fac;

            // Update (again) stmsta for next time step:
            stmsta = Math.max(0, stmsta - ((dry_wgt - old_dry_wgt) + meta)/assi);

            // Calculate std_len from from weight information:
            std_len = IBMFunction_NonEggStageBIOENGrowthRateDW.getL_fromW(dry_wgt, old_std_len, 2.928, -0.09, 0.047); // get parameters from R script

        }

        // Length and weight at t:
        //std_len += grSL*dtday; //mm dSL/dt = grSL
        // dry_wgt *= Math.exp(grDW*dtday); // mg dDW/dt = grDW*DW. The right equation is: Math.exp(grDW*dtday)
        //dry_wgt += gr_mg_fac;

        // Print light in grSL field, just to check the calculation (temporally):
        grSL = sum_ing;
        grDW = gr_mg_fac;

        // END OF BIOEN CHANGES ------------------------------------

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
     * Function to calculate movement rates.
     * 
     * @param dt - time step
     * @return 
     */
    public double[] calcUVW(double[] pos, double dt) {
        //compute vertical velocity
        double w = 0;
        if (typeVM==FDLStageParameters.FCN_VM_DVM_FixedDepthRanges) {
            //fcnVM instanceof wts.models.DisMELS.IBMFunctions.Movement.DielVerticalMigration_FixedDepthRanges
            //calculate the vertical swimming rate
                if(T<=0.0) T=0.01; 
                if (typeVV==FDLStageParameters.FCN_VV_FDL_VerticalSwimmingSpeed){
                    double TL = (std_len + 0.5169)/0.9315; //transform SL to TL
                    w = (Double) fcnVV.calculate(new Double[]{T,TL});//in mm/s
                    w = w/1000.0;//convert to m/s
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
            double[] res = (double[]) fcnVM.calculate(new double[]{dt,depth,td,w,90.833-ss[4]});
            w = res[0];
            attached = res[1]<0;
            if (attached) pos[2] = 0;//set individual on bottom
        }
        
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
        return new double[]{Math.signum(dt)*uv[0],Math.signum(dt)*uv[1],Math.signum(dt)*w};
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
        if (typeMort==FDLStageParameters.FCN_Mortality_ConstantMortalityRate){
            //fcnMortality instanceof ConstantMortalityRate
            mortalityRate = (Double)fcnMortality.calculate(null);
        } else 
        if (typeMort==FDLStageParameters.FCN_Mortality_InversePowerLawMortalityRate){
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
        atts.setValue(FDLStageAttributes.PROP_startTime,startTime);
        atts.setValue(FDLStageAttributes.PROP_time,time);
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
     * Updates attribute values defined for this class. 
     */
    @Override
    protected void updateAttributes() {
        super.updateAttributes();
        atts.setValue(FDLStageAttributes.PROP_attached,attached);
        atts.setValue(FDLStageAttributes.PROP_SL,std_len);
        atts.setValue(FDLStageAttributes.PROP_DW,dry_wgt);
        atts.setValue(FDLStageAttributes.PROP_stmsta,stmsta);
        atts.setValue(FDLStageAttributes.PROP_grSL,grSL);
        atts.setValue(FDLStageAttributes.PROP_grDW,grDW);
        atts.setValue(FDLStageAttributes.PROP_temperature,temperature);    
        atts.setValue(FDLStageAttributes.PROP_salinity,salinity);
        atts.setValue(FDLStageAttributes.PROP_rho,rho);
        atts.setValue(FDLStageAttributes.PROP_copepod,copepod);
        atts.setValue(FDLStageAttributes.PROP_euphausiid,euphausiid);
        atts.setValue(FDLStageAttributes.PROP_neocalanus,neocalanus);
    }

    /**
     * Updates local variables from the attributes.  
     */
    @Override
    protected void updateVariables() {
        super.updateVariables();
        attached     = atts.getValue(FDLStageAttributes.PROP_attached,attached);
        std_len      = atts.getValue(FDLStageAttributes.PROP_SL,std_len);
        dry_wgt      = atts.getValue(FDLStageAttributes.PROP_DW,dry_wgt);
        stmsta      = atts.getValue(FDLStageAttributes.PROP_stmsta,stmsta);
        grSL         = atts.getValue(FDLStageAttributes.PROP_grSL,grSL);
        grDW         = atts.getValue(FDLStageAttributes.PROP_grDW,grDW);
        temperature = atts.getValue(FDLStageAttributes.PROP_temperature,temperature);
        salinity    = atts.getValue(FDLStageAttributes.PROP_salinity,salinity);
        rho         = atts.getValue(FDLStageAttributes.PROP_rho,rho);
        copepod     = atts.getValue(FDLStageAttributes.PROP_copepod,copepod);
        euphausiid  = atts.getValue(FDLStageAttributes.PROP_euphausiid,euphausiid);
        neocalanus  = atts.getValue(FDLStageAttributes.PROP_neocalanus,neocalanus);
     }

}
