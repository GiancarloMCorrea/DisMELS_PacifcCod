/**
 * BenthicJuvStage.java
 *
 * Revisions: 
 * 20181011: 1. Corrected gL formula.
 *           2. Removed fcnDevelopment to match Parameters class.
 *           3. Removed fcnVV, fnVM because calculation for w is hard-wired to 0 in calcUVW.
 *           4. Added "attached" as new attribute (necessary with updated DisMELS).
 *           5. Removed "diam" since it's replaced by "length"
 * 20181015: 1. Removed minSettlementDepth, maxSettlementDepth parameters since they don't apply.
 *           2. Removed most of getMetamorphosedIndividuals() because there is no follow-on stage.
 * 20190722: 1. Removed fields associated with egg stage attributes "devStage" and "density"
 * 20210211: 1. Added DW, TL, and WW as attributes, with corresponding growth rates.
 *           2. Revised logic slightly for setAttributes, setInfo methods, updateAttributes, updateVariables
 *
 */

package sh.pcod.BenthicJuvStage;

import com.vividsolutions.jts.geom.Coordinate;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.List;
import org.openide.util.lookup.ServiceProvider;
import wts.models.DisMELS.IBMFunctions.Mortality.ConstantMortalityRate;
import wts.models.DisMELS.IBMFunctions.Mortality.InversePowerLawMortalityRate;
import wts.models.DisMELS.framework.*;
import wts.models.DisMELS.framework.IBMFunctions.IBMFunctionInterface;
import wts.roms.model.LagrangianParticle;
import sh.pcod.EpijuvStage.EpijuvStageAttributes;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageBIOENGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateSL;
import wts.models.DisMELS.IBMFunctions.HSMs.HSMFunction_Constant;
import wts.models.DisMELS.IBMFunctions.HSMs.HSMFunction_NetCDF;
import wts.models.DisMELS.IBMFunctions.HSMs.HSMFunction_NetCDF_InMemory;
import wts.models.utilities.CalendarIF;
import wts.roms.model.Interpolator3D;

/**
 * Class representing the P. cod benthic juvenile stage.
 * 
 * @author William Stockhausen
 * @author Sarah Hinckley
 * 
 */
@ServiceProvider(service=LifeStageInterface.class)
public class BenthicJuvStage extends AbstractLHS {
    
        //Static fields    
            //  Static fields new to this class
    /* flag to do debug operations */
    public static boolean debug = false;
    /* Class for attributes SH_NEW */
    public static final String attributesClass = 
            sh.pcod.BenthicJuvStage.BenthicJuvStageAttributes.class.getName();
    /* Class for parameters */
    public static final String parametersClass = 
            sh.pcod.BenthicJuvStage.BenthicJuvStageParameters.class.getName();
    /* Class for feature type for point positions */
    public static final String pointFTClass = 
            wts.models.DisMELS.framework.LHSPointFeatureType.class.getName();
    /* Classes for next LHS SH_NEW Add fdl*/
      public static final String[] nextLHSClasses = new String[]{ 
             sh.pcod.BenthicJuvStage.BenthicJuvStage.class.getName()};
    /* Classes for spawned LHS */
    public static final String[] spawnedLHSClasses = new String[]{};
    
    /* string identifying environmental field with copepod densities */
    private static final String FIELD_Cop = "Cop";
    /* string identifying environmental field with euphausiid densities */
    private static final String FIELD_Eup = "Eup";
    /* string identifying environmental field with neocalanus densities */
    private static final String FIELD_NCa = "NCa";
    /* string identifying environmental field with euphausiids shelf densities */
    private static final String EupS = "EupS";
    /* string identifying environmental field with microzooplankton densities */
    private static final String Mzl = "Mzl";
    /* string identifying environmental field with neocalanus shelf densities */
    private static final String NCaS = "NCaS";
    /* string identifying environmental field with large phytoplankton densities */
    private static final String PhL = "PhL";
    /* string identifying environmental field with small phytoplankton densities */
    private static final String PhS = "PhS";
    /* string identifying environmental field with u surface stress */
    private static final String Su = "Su";
    /* string identifying environmental field with v surface stress */
    private static final String Sv = "Sv";
    /* string identifying environmental field with u surface stress */
    private static final String pCO2 = "pCO2";
    /* string identifying environmental field with v surface stress */
    private static final String Mzs = "Mzs";
    
    //Instance fields
            //  Fields hiding ones from superclass
    /* life stage atrbutes object */
    protected BenthicJuvStageAttributes atts = null;
    /* life stage parameters object */
    protected BenthicJuvStageParameters params = null;
    
    //  Fields new to class
        //fields that reflect parameter values
    protected boolean isSuperIndividual;
    protected double  horizRWP;
    protected double  minStageDuration;
    protected double  maxStageDuration;
    protected double  stageTransRate;
    protected boolean useRandomTransitions;
    
        //fields that reflect (new) attribute values
    protected boolean attached = false;
    /** standard length (mm) */
    protected double std_len = 0;
    /** dry weight (mg) */
    protected double dry_wgt = 0;
    /** stomach state (units) */
    protected double stmsta = 0;
    /** stomach state (units) */
    protected double psurvival = 1;
    /** stomach state (units) */
    protected double mortfish = 0; // Here initial value of p_survival
    /** stomach state (units) */
    protected double mortinv = 0; // Here initial value of p_survival
    /** stomach state (units) */
    protected double avgRank = 0; // Here initial value of p_survival
    /** stomach state (units) */
    protected double avgSize = 0; // Here initial value of p_survival
    /** stomach state (units) */
    protected double pCO2val = 0;
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
     /** in situ neocalanoid (large copepods) density (mg/m^3, dry wt) */
    protected double neocalanus = 0;
     /** in situ neocalanoid (large copepods) density (mg/m^3, dry wt) */
    protected double neocalanusShelf = 0;
     /** in situ neocalanoid (large copepods) density (mg/m^3, dry wt) */
    protected double microzoo = 0;
    /** total length (mm) */
    protected double tot_len = 0;
    /** wet weight (mg) */
    protected double wet_wgt = 0;
    /** growth rate for total length (mm/d) */
    protected double grTL = 0;
    /** growth rate for wet weight (1/d) */
    protected double grWW = 0;
    /** habitat suitability index value */
    protected double hsi = 0;
        
            //other fields
    /** number of individuals transitioning to next stage */
    private double numTrans;  
    
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
    /** IBM function selected for growth in TL */
    private IBMFunctionInterface fcnGrTL = null; 
    /** IBM function selected for growth in WW */
    private IBMFunctionInterface fcnGrWW = null; 
    /** IBM function selected for HSM */
    private IBMFunctionInterface fcnHSI = null; 
    
    private int typeMort = 0;//integer indicating mortality function
    private int typeGrSL = 0;//integer indicating SL growth function
    private int typeGrDW = 0;//integer indicating DW growth function
    private int typeVM   = 0;//integer indicating vertical movement function
    private int typeVV   = 0;//integer indicating vertical velocity function
    private int typeGrTL = 0;//integer indicating TL growth function
    private int typeGrWW = 0;//integer indicating WW growth function
    private int typeHSI  = 0;//integer indicating HSI function

    private static final Logger logger = Logger.getLogger(BenthicJuvStage.class.getName());
    
    /**
     * Creates a new instance of BenthicJuvStage.  
     *  This constructor should be used ONLY to obtain
     *  the class names of the associated classes.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     */
    public BenthicJuvStage() {
        super("");
        super.atts = atts;
        super.params = params;
    }
    
    /**
     * Creates a new instance of BenthicJuvStage with the given typeName.
     * A new id number is calculated in the superclass and assigned to
     * the new instance's id, parentID, and origID. 
     * 
     * The attributes are vanilla.  Initial attribute values should be set,
     * then initialize() should be called to initialize all instance variables.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     * 
     * @param typeName
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     */
    public BenthicJuvStage(String typeName) 
                throws InstantiationException, IllegalAccessException {
        super(typeName);
        atts = new BenthicJuvStageAttributes(typeName);
        atts.setValue(LifeStageAttributesInterface.PROP_id,id);
        atts.setValue(LifeStageAttributesInterface.PROP_parentID,id);
        atts.setValue(LifeStageAttributesInterface.PROP_origID,id);
        params = (BenthicJuvStageParameters) LHS_Factory.createParameters(typeName);
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
     * 
     * 
     * @param strv - attributes as string array
     * @return - instance of LHS
     * @throws java.lang.InstantiationException
     * @throws java.lang.IllegalAccessException
     */
    @Override
    public BenthicJuvStage createInstance(String[] strv) 
                        throws InstantiationException, IllegalAccessException {
        LifeStageAttributesInterface theAtts = LHS_Factory.createAttributes(strv);
        BenthicJuvStage lhs = createInstance(theAtts);
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
    public BenthicJuvStage createInstance(LifeStageAttributesInterface theAtts)
                        throws InstantiationException, IllegalAccessException {
        BenthicJuvStage lhs = null;
        if (theAtts instanceof BenthicJuvStageAttributes) {
            lhs = new BenthicJuvStage(theAtts.getTypeName());
            long newID = lhs.id;//save id of new instance
            lhs.setAttributes(theAtts);
            if (lhs.atts.getID()==-1) {
                //constructing new individual, so reset id values to those of new
                lhs.id = newID;
                lhs.atts.setValue(BenthicJuvStageAttributes.PROP_id,newID);
            }
            newID = (Long) lhs.atts.getValue(BenthicJuvStageAttributes.PROP_parentID);
            if (newID==-1) {
                lhs.atts.setValue(BenthicJuvStageAttributes.PROP_parentID,newID);
            }
            newID = (Long) lhs.atts.getValue(BenthicJuvStageAttributes.PROP_origID);
            if (newID==-1) {
                lhs.atts.setValue(BenthicJuvStageAttributes.PROP_origID,newID);
            }
            lhs.initialize();//initialize instance variables
        }
        return lhs;
    }

    /**
     *  Returns the associated attributes.  
     * 
     * @return returns the attributes instance
     */
    @Override
    public BenthicJuvStageAttributes getAttributes() {
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
        aid = atts.getValue(BenthicJuvStageAttributes.PROP_id, id);
        if (aid==-1) {
            //change atts id to lhs id
            atts.setValue(BenthicJuvStageAttributes.PROP_id, id);
        } else {
            //change lhs id to atts id
            id = aid;
        }
        aid = atts.getValue(BenthicJuvStageAttributes.PROP_parentID, id);
        if (aid==-1) {
            atts.setValue(BenthicJuvStageAttributes.PROP_parentID, id);
        }
        aid = atts.getValue(BenthicJuvStageAttributes.PROP_origID, id);
        if (aid==-1) {
            atts.setValue(BenthicJuvStageAttributes.PROP_origID, id);
        }
        initialize();//initialize instance variables
    }

    /**
     * Sets the attributes for the instance by copying values from the input.
     * This does NOT change the typeName of the LHS instance (or the associated 
     * LHSAttributes instance) on which the method is called.
     * Note that ALL attributes are copied, so id, parentID, and origID are copied
     * as well. 
     * <pre>
     *  Side effects:
     *      1. updateVariables() is called to update instance variables.
     *      2. Instance field "id" is also updated.
     * </pre>
     * @param newAtts - should be instance of SimplePelagicLHSAttributes
     */
    @Override
    public void setAttributes(LifeStageAttributesInterface newAtts) {
        if (newAtts instanceof BenthicJuvStageAttributes) {
            BenthicJuvStageAttributes oldAtts = (BenthicJuvStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));}
        else if(newAtts instanceof EpijuvStageAttributes) {
            EpijuvStageAttributes oldAtts = (EpijuvStageAttributes) newAtts;
            //EpijuvStageAttributes and BenthicJuvStageAttributes have identical attribute keys 
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));
            //no need to set attributes NOT in EpijuvStageAttributes
        } else {
            //TODO: should throw an error here
            logger.info("setAttributes(): no match for attributes type:"+newAtts.toString());
        }
        id = atts.getValue(BenthicJuvStageAttributes.PROP_id, id);
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
        atts.setValue(BenthicJuvStageAttributes.PROP_ageInStage, 0.0);//reset age in stage
        atts.setValue(BenthicJuvStageAttributes.PROP_active,true);    //set active to true
        atts.setValue(BenthicJuvStageAttributes.PROP_alive,true);     //set alive to true
        id = atts.getID(); //reset id for current LHS to one from old LHS

        //update local variables to capture changes made here
        updateVariables();
    }
    
    /**
     *  Sets the associated attributes object. Use this after creating an LHS instance
     * as an "output" from another LHS that is functioning as a super individual.
     * 
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
        atts.setValue(BenthicJuvStageAttributes.PROP_id,idc);//reset id to one for current LHS
        atts.setValue(BenthicJuvStageAttributes.PROP_parentID,
                      oldAtts.getValue(LifeStageAttributesInterface.PROP_id));//copy old id to parentID
        atts.setValue(BenthicJuvStageAttributes.PROP_number, numTrans);//set number to numTrans
        atts.setValue(BenthicJuvStageAttributes.PROP_ageInStage, 0.0); //reset age in stage
        atts.setValue(BenthicJuvStageAttributes.PROP_active,true);     //set active to true
        atts.setValue(BenthicJuvStageAttributes.PROP_alive,true);      //set alive to true
            
        //update local variables to capture changes made here
        updateVariables();
    }

    /**
     *  Returns the associated parameters.  
     */
    @Override
    public BenthicJuvStageParameters getParameters() {
        return params;
    }

    /**
     * Sets the parameters for the instance to a cloned version of the input.
     * @param newParams - should be instance of BenthicJuvStageParameters
     */
    @Override
    public void setParameters(LifeStageParametersInterface newParams) {
        if (newParams instanceof BenthicJuvStageParameters) {
            params = (BenthicJuvStageParameters) newParams;
            super.params = params;
            setParameterValues();
            fcnMortality = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_Mortality);
            fcnGrSL = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_GrowthSL);
            fcnGrDW = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_GrowthDW);
            fcnGrTL = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_GrowthTL);
            fcnGrWW = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_GrowthWW);
            fcnHSI  = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_HSM);
            
            if (fcnMortality instanceof ConstantMortalityRate)
                typeMort = BenthicJuvStageParameters.FCN_Mortality_ConstantMortalityRate;
            else if (fcnMortality instanceof InversePowerLawMortalityRate)
                typeMort = BenthicJuvStageParameters.FCN_Mortality_InversePowerLawMortalityRate;
            
            if (fcnGrSL instanceof IBMFunction_NonEggStageSTDGrowthRateSL) 
                typeGrSL = BenthicJuvStageParameters.FCN_GrSL_NonEggStageSTDGrowthRate;
            
            if (fcnGrDW instanceof IBMFunction_NonEggStageSTDGrowthRateDW) 
                typeGrDW = BenthicJuvStageParameters.FCN_GrDW_NonEggStageSTDGrowthRate;
            else if (fcnGrDW instanceof IBMFunction_NonEggStageBIOENGrowthRateDW)
                typeGrDW = BenthicJuvStageParameters.FCN_GrDW_NonEggStageBIOENGrowthRate;
            
            if (fcnGrTL instanceof IBMFunction_BenthicJuv_GrowthRateTL)    
                typeGrTL = BenthicJuvStageParameters.FCN_GrTL_BenthicJuv_GrowthRate;
            
            if (fcnGrWW instanceof IBMFunction_BenthicJuv_GrowthRateWW)    
                typeGrWW = BenthicJuvStageParameters.FCN_GrWW_BenthicJuv_GrowthRate;
            
            if (fcnHSI instanceof HSMFunction_Constant)        
                typeHSI = BenthicJuvStageParameters.FCN_HSM_Constant;
            else if (fcnHSI instanceof HSMFunction_NetCDF)          
                typeHSI = BenthicJuvStageParameters.FCN_HSM_NetCDF;
            else if (fcnHSI instanceof HSMFunction_NetCDF_InMemory) 
                typeHSI = BenthicJuvStageParameters.FCN_HSM_NetCDF_InMemory;
        } else {
            //TODO: throw some error
        }
    }
    
    /*
     * Copy the values from the params map to the param variables.
     */
    private void setParameterValues() {
        isSuperIndividual = 
                params.getValue(BenthicJuvStageParameters.PARAM_isSuperIndividual,isSuperIndividual);
        horizRWP = 
                params.getValue(BenthicJuvStageParameters.PARAM_horizRWP,horizRWP);
        minStageDuration = 
                params.getValue(BenthicJuvStageParameters.PARAM_minStageDuration,minStageDuration);
        maxStageDuration = 
                params.getValue(BenthicJuvStageParameters.PARAM_maxStageDuration,maxStageDuration);
        useRandomTransitions = 
                params.getValue(BenthicJuvStageParameters.PARAM_useRandomTransitions,true);
    }
    
    /**
     *  Provides a copy of the object.  The attributes and parameters
     *  are cloned in the process, so the clone is independent of the
     *  original.
     */
    @Override
    public Object clone() {
        BenthicJuvStage clone = null;
        try {
            clone = (BenthicJuvStage) super.clone();
            clone.setAttributes(atts);  //this clones atts
            clone.updateVariables();    //this sets the variables in the clone to the attribute values
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
     * No "next" life stage for BenthicJuv individuals, 
     * so no metamorphosed individuals.
     * 
     * @param dt - time step in seconds
     * @return 
     */
    @Override
    public List<LifeStageInterface> getMetamorphosedIndividuals(double dt) {
        output.clear();
        return output;
    
    }
    
    /**
     * This would ordinarily be called from getMetamorphosedIndividuals(), but
     * BenthicJuvs do not have a "next" life stage.
     * 
     * @return 
     */
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
        hType      = atts.getValue(BenthicJuvStageAttributes.PROP_horizType,hType);
        vType      = atts.getValue(BenthicJuvStageAttributes.PROP_vertType,vType);
        xPos       = atts.getValue(BenthicJuvStageAttributes.PROP_horizPos1,xPos);
        yPos       = atts.getValue(BenthicJuvStageAttributes.PROP_horizPos2,yPos);
        zPos       = atts.getValue(BenthicJuvStageAttributes.PROP_vertPos,zPos);
        time       = startTime;
        numTrans   = 0.0; //set numTrans to zero
        if (debug) logger.info(hType+cc+vType+cc+startTime+cc+xPos+cc+yPos+cc+zPos);
        if (i3d!=null) {
            double[] IJ = new double[] {xPos,yPos};
            if (hType==Types.HORIZ_XY) {
                IJ = i3d.getGrid().computeIJfromXY(xPos,yPos);
            } else if (hType==Types.HORIZ_LL) {
//                if (xPos<0) xPos=xPos+360;
                IJ = i3d.getGrid().computeIJfromLL(yPos,xPos);
            }
            double z = i3d.interpolateBathymetricDepth(IJ);
            if (debug) logger.info("Bathymetric depth = "+z);
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
            atts.setValue(BenthicJuvStageAttributes.PROP_horizType,Types.HORIZ_LL);
            atts.setValue(BenthicJuvStageAttributes.PROP_vertType,Types.VERT_H);
            
            //interpolate initial position and environmental variables
            double[] pos = lp.getIJK();
            updatePosition(pos);
            updateEnvVars(pos);
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
        //SH_NEW
        double T = i3d.interpolateTemperature(pos);
        if(T<=0.0) T=0.01; 
        
        //SH-Prey Stuff  
        copepod    = i3d.interpolateValue(pos,FIELD_Cop,Interpolator3D.INTERP_VAL);
        euphausiid = i3d.interpolateValue(pos,FIELD_Eup,Interpolator3D.INTERP_VAL);
        neocalanus = i3d.interpolateValue(pos,FIELD_NCa,Interpolator3D.INTERP_VAL);
        neocalanusShelf  = i3d.interpolateValue(pos,NCaS,Interpolator3D.INTERP_VAL);
        double euphausiids_shelf  = i3d.interpolateValue(pos,EupS,Interpolator3D.INTERP_VAL);
        double microzoo_large  = i3d.interpolateValue(pos,Mzl,Interpolator3D.INTERP_VAL);
        double phytoL = i3d.interpolateValue(pos,PhL,Interpolator3D.INTERP_VAL);
        double phytoS = i3d.interpolateValue(pos,PhS,Interpolator3D.INTERP_VAL);
        double tauX = i3d.interpolateValue(pos2d,Su,Interpolator3D.INTERP_VAL); // 3D interpolator but should use 2D internally
        double tauY = i3d.interpolateValue(pos2d,Sv,Interpolator3D.INTERP_VAL); // 3D interpolator but should use 2D internally
        double chlorophyll = (phytoL/25) + (phytoS/65); // calculate chlorophyll (mg/m^-3) 
        // values 25 and 65 based on Kearney et al 2018 Table A4
        double microzoo_small  = i3d.interpolateValue(pos,Mzs,Interpolator3D.INTERP_VAL);
        microzoo = microzoo_small + microzoo_large; // total microzooplankton
        pCO2val  = i3d.interpolateValue(pos,pCO2,Interpolator3D.INTERP_VAL);
        //calculate pco2 factor:
        double facCO2 = 0;
        if(pCO2val > 600) { // only for values larger than 600
            facCO2 = IBMFunction_NonEggStageBIOENGrowthRateDW.calcCO2(pCO2val);
        } 

        time += dt;
        double dtday = dt/86400;//time step in days

        // Create objects for output of BIOEN:
        double gr_mg_fac = 0; // factor to avoid dry_wgt in IF 
        double meta = 0;
        double costRateOfMetabolism = 0.5; // check this number
        double old_dry_wgt = dry_wgt; // save previous dry_wgt
        double old_std_len = std_len;
        // Light (begin):
        // create object for light calculation:
        double[] eb2 = new double[2]; // K parameter + second part of Eb equation
        double eb = 0; // create Eb object
        CalendarIF cal2 = null;
        double[] ltemp = new double[3];
        double[] ltemp2 = new double[2];
        cal2 = GlobalInfo.getInstance().getCalendar(); // to calculate julian day
        ltemp = IBMFunction_NonEggStageBIOENGrowthRateDW.calcLightQSW(lat,cal2.getYearDay()); // see line 713 in ibm.py
        double maxLight = ltemp[0]/0.217; // see line 714 in ibm.py
        ltemp2 = IBMFunction_NonEggStageBIOENGrowthRateDW.calcLightSurlig(lat,cal2.getYearDay(), maxLight); // see line 715 in ibm.py
        eb2 = IBMFunction_NonEggStageBIOENGrowthRateDW.calcLight(chlorophyll, depth, bathym); // K parameter and second part of Eb equation
        eb = 0.42*ltemp2[1]*eb2[1]; // see line 727 in ibm.py. This is Eb. 0.42 as in Kearney et al 2020 Eq A14
        // Light (end):

        //calculate growth in length, weight
        if(T<=0.0) T=0.01; 
        if (typeGrSL==BenthicJuvStageParameters.FCN_GrSL_NonEggStageSTDGrowthRate) {
            grSL = (Double) fcnGrSL.calculate(new Double[]{T,std_len});
            std_len += grSL*dtday;
        }
        if (typeGrDW==BenthicJuvStageParameters.FCN_GrDW_NonEggStageSTDGrowthRate) {
            grDW = (Double) fcnGrDW.calculate(new Double[]{T,dry_wgt});
            dry_wgt *= Math.exp(grDW * dtday);
        }
        if (typeGrTL==BenthicJuvStageParameters.FCN_GrTL_BenthicJuv_GrowthRate) {
            grTL = (Double) fcnGrTL.calculate(T);
            tot_len += grTL*dtday;
        }

        if (typeGrDW==BenthicJuvStageParameters.FCN_GrDW_NonEggStageBIOENGrowthRate){

            double r = ((0.454 + 1.610*T - 0.069*T*T)*Math.exp(-6.725*old_dry_wgt))/100;// units: 1/day. This is similar to 'g' (1/day) in TROND. 
            grDW = old_dry_wgt*(Math.exp(r*dtday) - 1); // Same as TROND

            // Begin function:
            meta = dtday*2.38e-7*Math.exp(0.088*T)*Math.pow(old_dry_wgt,0.9)*(1 + facCO2*0.1); // as in Kristiansen et al 2007. Units: mg/day (without dt). HERE I CHANGED dt FOR dtday 
            // dtday makes more sense 

            if(eb > 0.001) {
                if(std_len > 5.5){
                    meta *= 2.5;
                } else {
                    meta *= 1.4;
                }
            } 

            double activityCost = 0.5*meta*costRateOfMetabolism; // TODO: (diffZ/maxDiffZ) = 0.5, but this should change based on vertical movement
            gr_mg_fac = grDW - activityCost; // Here grDW is as gr_mg in TROND

            // New weight in mg:
            dry_wgt += gr_mg_fac;

            // Calculate std_len from from weight information:
            std_len = IBMFunction_NonEggStageBIOENGrowthRateDW.getL_fromW(dry_wgt, old_std_len); // get parameters from R script
            grSL = std_len - old_std_len;

            tot_len = std_len*1.025; // Just multiply by a factor, TODO: discuss this later

        }

        if (typeGrWW==BenthicJuvStageParameters.FCN_GrWW_BenthicJuv_GrowthRate)
            grWW = (Double) fcnGrWW.calculate(T);
        wet_wgt *= Math.exp(grWW * dtday);
        
        // Survival rate (begin):
        psurvival = 0;
        mortfish = 0;
        mortinv = 0;
        avgRank = 0;
        avgSize = 0;
        // Survival rate (end):

        updatePosition(pos);
        updateEnvVars(pos);
        updateNum(dt);
        updateAge(dt);
        
        //check for exiting grid (no real need for this until BenthicJuvs move)
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
        if (typeMort == BenthicJuvStageParameters.FCN_Mortality_ConstantMortalityRate){
            mortalityRate = (Double)fcnMortality.calculate(null);
        } else 
        if (typeMort == BenthicJuvStageParameters.FCN_Mortality_InversePowerLawMortalityRate){
          mortalityRate = (Double)fcnMortality.calculate(std_len);//using standard length as covariate for mortality
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
    }
    
    /**
     * Updates location-associated information given the current location.
     * <pre>
     * Updates the following position-related variables:
     *   1. bathymetry (bathym)
     *   2. depth
     *   3. latitude (lat)
     *   4. longitude (lon)
     *   5. gridCellID
     *   6. the track
     * </pre>
     * @param pos - double[] giving position in ROMS {xi, eta, K} grid coordinates
     */
    private void updatePosition(double[] pos) {
        bathym     =  i3d.interpolateBathymetricDepth(pos);
        depth      = -i3d.calcZfromK(pos[0],pos[1],pos[2]);
        lat        =  i3d.interpolateLat(pos);
        lon        =  i3d.interpolateLon(pos);
        gridCellID = ""+Math.round(pos[0])+"_"+Math.round(pos[1]);
        updateTrack();
    }
    
    /**
     * Interpolates the environmental variables to the current location.
     * <pre>
     * Updates the following environmental variables:
     *   1. temperature
     *   2. salinity
     *   3. seawater density (rho, if included in the physical environment)
     *   4. copepod density (if included in the physical environment)
     *   5. euphausiid density (if included in the physical environment)
     *   6. neocalanus density (if included in the physical environment)
     *   7. hsi
     * </pre>
     * @param pos - double[] giving position in ROMS {xi, eta, K} grid coordinates
     */
    private void updateEnvVars(double[] pos) {
        temperature = i3d.interpolateTemperature(pos);
        salinity    = i3d.interpolateSalinity(pos);
        
        if (i3d.getPhysicalEnvironment().getField("rho")!=null) 
            rho  = i3d.interpolateValue(pos,"rho");
        else rho = 0.0;
        
        if (i3d.getPhysicalEnvironment().getField(FIELD_Cop)!=null) 
            copepod    = i3d.interpolateValue(pos,FIELD_Cop,Interpolator3D.INTERP_VAL);
        
        if (i3d.getPhysicalEnvironment().getField(FIELD_Eup)!=null) 
            euphausiid = i3d.interpolateValue(pos,FIELD_Eup,Interpolator3D.INTERP_VAL);
        
        if (i3d.getPhysicalEnvironment().getField(FIELD_NCa)!=null) 
            neocalanus = i3d.interpolateValue(pos,FIELD_NCa,Interpolator3D.INTERP_VAL);
        
        switch (typeHSI) {
            case BenthicJuvStageParameters.FCN_HSM_Constant:
                hsi = (Double)fcnHSI.calculate(null);//constant value
                break;
            case BenthicJuvStageParameters.FCN_HSM_NetCDF:
                {
                    double[] posLL = new double[]{lon,lat};
                    hsi = (Double)fcnHSI.calculate(posLL);
                    break;
                }
            case BenthicJuvStageParameters.FCN_HSM_NetCDF_InMemory:
                {
                    double[] posLL = new double[]{lon,lat};
                    hsi = (Double)fcnHSI.calculate(posLL);
                    break;
                }
            default:
                break;
        }
    }

    @Override
    public double getStartTime() {
        return startTime;
    }

    @Override
    public void setStartTime(double newTime) {
        startTime = newTime;
        time      = startTime;
        atts.setValue(BenthicJuvStageAttributes.PROP_startTime,startTime);
        atts.setValue(BenthicJuvStageAttributes.PROP_time,     time);
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
        atts.setValue(BenthicJuvStageAttributes.PROP_attached,   attached);
        atts.setValue(BenthicJuvStageAttributes.PROP_SL,         std_len);
        atts.setValue(BenthicJuvStageAttributes.PROP_DW,         dry_wgt);
        atts.setValue(BenthicJuvStageAttributes.PROP_stmsta,     stmsta);
        atts.setValue(BenthicJuvStageAttributes.PROP_psurvival,  psurvival);
        atts.setValue(BenthicJuvStageAttributes.PROP_mortfish,mortfish);
        atts.setValue(BenthicJuvStageAttributes.PROP_mortinv,mortinv);
        atts.setValue(BenthicJuvStageAttributes.PROP_avgRank,avgRank);
        atts.setValue(BenthicJuvStageAttributes.PROP_avgSize,avgSize);
        atts.setValue(BenthicJuvStageAttributes.PROP_pCO2val,  pCO2val);
        atts.setValue(BenthicJuvStageAttributes.PROP_grSL,       grSL);
        atts.setValue(BenthicJuvStageAttributes.PROP_grDW,       grDW);
        atts.setValue(BenthicJuvStageAttributes.PROP_temperature,temperature);    
        atts.setValue(BenthicJuvStageAttributes.PROP_salinity,   salinity);
        atts.setValue(BenthicJuvStageAttributes.PROP_rho,        rho);
        atts.setValue(BenthicJuvStageAttributes.PROP_copepod,    copepod);
        atts.setValue(BenthicJuvStageAttributes.PROP_neocalanus, neocalanus);
        atts.setValue(BenthicJuvStageAttributes.PROP_euphausiid, euphausiid);
        atts.setValue(BenthicJuvStageAttributes.PROP_neocalanusShelf, neocalanusShelf);
        atts.setValue(BenthicJuvStageAttributes.PROP_microzoo, microzoo);
        atts.setValue(BenthicJuvStageAttributes.PROP_TL,         tot_len);
        atts.setValue(BenthicJuvStageAttributes.PROP_WW,         wet_wgt);
        atts.setValue(BenthicJuvStageAttributes.PROP_grTL,       grTL);
        atts.setValue(BenthicJuvStageAttributes.PROP_grWW,       grWW);
        atts.setValue(BenthicJuvStageAttributes.PROP_hsi,        hsi);

    }

    /**
     * Updates local variables from the attributes.  
     */
    @Override
    protected void updateVariables() {
        super.updateVariables();
        attached    = atts.getValue(BenthicJuvStageAttributes.PROP_attached,    attached);
        std_len     = atts.getValue(BenthicJuvStageAttributes.PROP_SL,          std_len);
        dry_wgt     = atts.getValue(BenthicJuvStageAttributes.PROP_DW,          dry_wgt);
        stmsta      = atts.getValue(BenthicJuvStageAttributes.PROP_stmsta,      stmsta);
        psurvival    = atts.getValue(BenthicJuvStageAttributes.PROP_psurvival,  psurvival);
        mortfish      = atts.getValue(BenthicJuvStageAttributes.PROP_mortfish,mortfish);
        mortinv      = atts.getValue(BenthicJuvStageAttributes.PROP_mortinv,mortinv);
        avgRank      = atts.getValue(BenthicJuvStageAttributes.PROP_avgRank,avgRank);
        avgSize      = atts.getValue(BenthicJuvStageAttributes.PROP_avgSize,avgSize);
        pCO2val    = atts.getValue(BenthicJuvStageAttributes.PROP_pCO2val,  pCO2val);
        grSL        = atts.getValue(BenthicJuvStageAttributes.PROP_grSL,        grSL);
        grDW        = atts.getValue(BenthicJuvStageAttributes.PROP_grDW,        grDW);
        temperature = atts.getValue(BenthicJuvStageAttributes.PROP_temperature, temperature);
        salinity    = atts.getValue(BenthicJuvStageAttributes.PROP_salinity,    salinity);
        rho         = atts.getValue(BenthicJuvStageAttributes.PROP_rho,         rho);
        copepod     = atts.getValue(BenthicJuvStageAttributes.PROP_copepod,     copepod);
        neocalanus  = atts.getValue(BenthicJuvStageAttributes.PROP_neocalanus,  neocalanus);
        euphausiid  = atts.getValue(BenthicJuvStageAttributes.PROP_euphausiid,  euphausiid);
        neocalanusShelf  = atts.getValue(BenthicJuvStageAttributes.PROP_neocalanusShelf,  neocalanusShelf);
        microzoo  = atts.getValue(BenthicJuvStageAttributes.PROP_microzoo,  microzoo);
        tot_len     = atts.getValue(BenthicJuvStageAttributes.PROP_TL,          tot_len);
        wet_wgt     = atts.getValue(BenthicJuvStageAttributes.PROP_WW,          wet_wgt);
        grTL        = atts.getValue(BenthicJuvStageAttributes.PROP_grTL,        grTL);
        grWW        = atts.getValue(BenthicJuvStageAttributes.PROP_grWW,        grWW);
        hsi         = atts.getValue(BenthicJuvStageAttributes.PROP_hsi,         hsi);
     }
}
