/*
 * Postflexion Pcod Larval Stage
 * Devstage = 3
 * Created on January 24, 2006, 11:33 AM
 *
 * To change this template, choose Tools | Template Manager
 * and open the template in the editor.
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
import wts.models.utilities.DateTimeFunctions;
import wts.models.utilities.ModelCalendar;
import wts.roms.model.LagrangianParticle;
import sh.pcod.FDLpfStage.FDLpfStageAttributes;
import sh.pcod.EpijuvStage.EpijuvStageAttributes;
import wts.roms.model.Interpolator3D;

/**
 *
 * @author William Stockhausen
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
    protected double  minStageSize;
    protected double minSettlementDepth;
    protected double maxSettlementDepth;
    protected double  stageTransRate;
    protected boolean useRandomTransitions;
    
        //fields that reflect (new) attribute values
    /** development stage,0=egg,1=ysl,2=fdl,3=epijuvFDLpf,4=Epijuv, 5=BenthicJuv */
    protected double devStage;
    /** egg diameter (mm) */
    protected double diam = 0;
    /** density of egg [kg/m^3]) */
    protected double density;
    /** in situ temperature (deg C) */
    protected double temperature = 0;
    /** in situ salinity */
    protected double salinity = 0;
    /** in situ water density */
    protected double rho = 0;
     protected double copepod = 0;    /** in situ small copepods */
     /** in situ euphausiid density mg/m^3, dry wt)) */
    protected double euphausiid = 0;    /** in situ euphausiids */
     /** in situ neocalanoid density mg/m^3, dry wt)) */
    protected double neocalanus = 0;    /** in situ large copepods */      
    /**growth in Length mm/d */
    /**SH_NEW*/
    protected double gL = 0;
        /**FDL Length variable (mm) */
    protected double length = 0;
    /**FDL maximum size = random between 25-35.  Stays the same at each time step*/
    protected double maxlength;
    //FDL Size at flexion
    //protected double flexion;
    double T0, T1, T;

            //other fields
    /** number of individuals transitioning to next stage */
    private double numTrans;  
    /** total depth (m) at individual's position */
    private double totalDepth;
    
    /** IBM function selected for development */
    private IBMFunctionInterface fcnDevelopment = null; 
    /** IBM function selected for mortality */
    private IBMFunctionInterface fcnMortality = null; 
    /** IBM function selected for vertical movement */
    private IBMFunctionInterface fcnVM = null; 
    /** IBM function selected for vertical velocity */
    private IBMFunctionInterface fcnVV = null; 
    
    private static final Logger logger = Logger.getLogger(BenthicJuvStage.class.getName());
    
    /**
     * Creates a new instance of GenericLHS.  
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
     * Creates a new instance of SimplePelagicLHS with the given typeName.
     * A new id number is calculated in the superclass and assigned to
     * the new instance's id, parentID, and origID. 
     * 
     * The attributes are vanilla.  Initial attribute values should be set,
     * then initialize() should be called to initialize all instance variables.
     * DO NOT DELETE THIS CONSTRUCTOR!!
     */
    public BenthicJuvStage(String typeName) 
                throws InstantiationException, IllegalAccessException {
        super(typeName);
        atts = new BenthicJuvStageAttributes(typeName);
        atts.setValue(BenthicJuvStageAttributes.PROP_id,id);
        atts.setValue(BenthicJuvStageAttributes.PROP_parentID,id);
        atts.setValue(BenthicJuvStageAttributes.PROP_origID,id);
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
     * /
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
        }
        lhs.initialize();//initialize instance variables
        return lhs;
    }

    /**
     *  Returns the associated attributes.  
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
     *  Side effects:
     *      updateVariables() is called to update instance variables.
     *      Instance field "id" is also updated.
     * @param newAtts - should be instance of SimplePelagicLHSAttributes
     */
    @Override
    public void setAttributes(LifeStageAttributesInterface newAtts) {
        if (newAtts instanceof BenthicJuvStageAttributes) {
            BenthicJuvStageAttributes oldAtts = (BenthicJuvStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));}
        else if(newAtts instanceof EpijuvStageAttributes) {
            EpijuvStageAttributes oldAtts = (EpijuvStageAttributes) newAtts;
            for (String key: atts.getKeys()) atts.setValue(key,oldAtts.getValue(key));
           //SH_NEW
            // atts.setValue(atts.PROP_length,oldAtts.getValue(EggStageAttributes.PROP_diameter, 1.0));
            atts.setValue(BenthicJuvStageAttributes.PROP_length,oldAtts.getValue(EpijuvStageAttributes.PROP_length, length));

        
        } else {
            
            
            
            //TODO: should throw an error here
            logger.info("setAttributes(): no match for attributes type:"+newAtts.toString());
        }
        id = atts.getValue(BenthicJuvStageAttributes.PROP_id, id);
        updateVariables();
    }
    
    /**
     *  Sets the associated attributes object. Use this after creating an LHS instance
     * as an "output" from another LHS that is functioning as an ordinary individual.
     */
    @Override
    public void setInfoFromIndividual(LifeStageInterface oldLHS){
        /** 
         * Since this is a single individual making a transition, we need to:
         *  1) copy the attributes from the old LHS (id's should remain as for old LHS)
         *  2) set age in stage = 0
         *  3) set active and alive to true
         *  5) copy the Lagrangian Particle from the old LHS
         *  6) start a new track from the current position for the oldLHS
         *  7) update local variables
         */
        LifeStageAttributesInterface oldAtts = oldLHS.getAttributes();            
        setAttributes(oldAtts);
        
        //reset some attributes
        atts.setValue(BenthicJuvStageAttributes.PROP_ageInStage, 0.0);//reset age in stage
        atts.setValue(BenthicJuvStageAttributes.PROP_active,true);    //set active to true
        atts.setValue(BenthicJuvStageAttributes.PROP_alive,true);     //set alive to true
        id = atts.getID(); //reset id for current LHS to one from old LHS

        //copy LagrangianParticle information
        this.setLagrangianParticle(oldLHS.getLagrangianParticle());
        //start track at last position of oldLHS track
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_PROJECTED),COORDINATE_TYPE_PROJECTED);
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_GEOGRAPHIC),COORDINATE_TYPE_GEOGRAPHIC);
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
         *          1) copy most attribute values from old stage
         *          2) make sure id for this LHS is retained, not changed
         *          3) assign old LHS id to this LHS as parentID
         *          4) copy old LHS origID to this LHS origID
         *          5) set number in this LHS to numTrans
         *          6) reset age in stage to 0
         *          7) set active and alive to true
         *          9) copy the Lagrangian Particle from the old LHS
         *         10) start a new track from the current position for the oldLHS
         *         11) update local variables to match attributes
         */
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
            
        //copy LagrangianParticle information
        this.setLagrangianParticle(oldLHS.getLagrangianParticle());
        //start track at last position of oldLHS track
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_PROJECTED),COORDINATE_TYPE_PROJECTED);
        this.startTrack(oldLHS.getLastPosition(COORDINATE_TYPE_GEOGRAPHIC),COORDINATE_TYPE_GEOGRAPHIC);
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
            fcnDevelopment = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_Development);
            fcnMortality = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_Mortality);
            fcnVM = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_VerticalMovement);
            fcnVV = params.getSelectedIBMFunctionForCategory(BenthicJuvStageParameters.FCAT_VerticalVelocity);
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
        minSettlementDepth = 
                params.getValue(BenthicJuvStageParameters.PARAM_minSettlementDepth,minSettlementDepth);
        maxSettlementDepth = 
                params.getValue(BenthicJuvStageParameters.PARAM_maxSettlementDepth,maxSettlementDepth);
        
        
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
            clone.setAttributes(atts);//this clones atts
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
        output.clear();
        LifeStageInterface nLHS = null;
        //if total depth is appropriate for settlement and 
        //indiv is near the bottom, then settle and transform to next stage.
        if ((totalDepth>=minSettlementDepth)&&
                (totalDepth<=maxSettlementDepth)&&
                (depth>(totalDepth-5))) {
                    devStage = 5;
                    //density = totalDepth;
             if ((numTrans>0)||!isSuperIndividual){
                nLHS = createNextLHS();
                if (nLHS!=null) output.add(nLHS);
                }
            }    
        return output;
    
    }
    
    
   /* @Override
    public List<LifeStageInterface> getMetamorphosedIndividuals(double dt) {
        double dtp = 0.25*(dt/86400);//use 1/4 timestep (converted from sec to d)
        output.clear();
        LifeStageInterface nLHS;
        
        //if(length >= flexion) {
          if(length >= 35.0) {
           devStage = 4;   
           if ((numTrans>0)||!isSuperIndividual){
                nLHS = createNextLHS();
                if (nLHS!=null) output.add(nLHS);
            }
        }

        return output;
    }*/
    
    private LifeStageInterface createNextLHS() {
        LifeStageInterface nLHS = null;
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
                nLHS = LHS_Factory.createNextLHSFromSuperIndividual(typeName,this,numTrans);
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
                nLHS = LHS_Factory.createNextLHSFromIndividual(typeName,this);
                alive  = false; //allow only 1 transition, so kill this stage
                active = false; //set stage inactive, also
            }
        } catch (IllegalAccessException | InstantiationException ex) {
            ex.printStackTrace();
        }
        return nLHS;
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
        hType      = atts.getValue(BenthicJuvStageAttributes.PROP_horizType,hType);
        vType      = atts.getValue(BenthicJuvStageAttributes.PROP_vertType,vType);
        xPos       = atts.getValue(BenthicJuvStageAttributes.PROP_horizPos1,xPos);
        yPos       = atts.getValue(BenthicJuvStageAttributes.PROP_horizPos2,yPos);
        zPos       = atts.getValue(BenthicJuvStageAttributes.PROP_vertPos,zPos);
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
            atts.setValue(BenthicJuvStageAttributes.PROP_horizType,Types.HORIZ_LL);
            atts.setValue(BenthicJuvStageAttributes.PROP_vertType,Types.VERT_H);
            //interpolate initial position and environmental variables
            double[] pos = lp.getIJK();
            //SH_NEW  9_14
            //updatePosition(pos);
            interpolateEnvVars(pos);
            updateAttributes(); 
        }
    }
    
    @Override
    public void step(double dt) throws ArrayIndexOutOfBoundsException {
        //WTS_NEW 2012-07-26:{
        double[] pos = lp.getIJK();
        //SH_NEW
        T0 = i3d.interpolateTemperature(pos);
        T1 = i3d.interpolateTemperature(pos);
        T = 0.5 * (T0 + T1);
        
             //SH-Prey Stuff  
        String Cop = "Cop";
        copepod = i3d.interpolateValue(pos,Cop,Interpolator3D.INTERP_VAL);
        String Eup = "Eup";
        euphausiid = i3d.interpolateValue(pos,Eup,Interpolator3D.INTERP_VAL);
        String NCa = "NCa";
        neocalanus = i3d.interpolateValue(pos,NCa,Interpolator3D.INTERP_VAL);
      
        
        
            
        if (attached){
            lp.setIJK(pos[0], pos[1], pos[2]);
        } else {
            double[] uvw = calcUVW(pos,dt);//this also sets "attached" and may change pos[2] to 0

            //}:WTS_NEW 2012-07-26
            //do lagrangian particle tracking
            lp.setU(uvw[0],lp.getN());
            lp.setV(uvw[1],lp.getN());
            lp.setW(uvw[2],lp.getN());
            //now do predictor step
            lp.doPredictorStep();
            //assume same daytime status, but recalc depth and revise W 
            pos = lp.getPredictedIJK();
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
        time = time+dt;
        
        
        //TODO: need to update devStage, diam, density, number
        devStage = 5;
        
        //SH_NEW:{
        double dtday = dt/86400;        //dt=biolmodel time step. At 72/day, dt(sec)= 1200; dtday=0.014
        //flexion = trian(10.0,13.5,17.0); // ??? What should I use here
        //flexion = 10.0;
        //Growth in length
        if(T<=0.0) T=0.01; 
        gL = (0.034 + (0.043 * T) - (0.0008 * T * T));//Hurst et al 2010, postflexion eq, mm per day
        length = length + (gL*dtday);
        
        
        //Need Mortality
        
        //No DVM until after flexion - next stage
        //Need Swimspeed
        //Preflexion larval depth: 0-40m

        diam = length;
        //density = length;

        //if (length >= flexion) devStage=3;        
        //System.out.print("stagedurYSL = " + stagedur);
        // Use diam variable to test my new algorithms as it is output...
        //}:SH_NEW
        
        updateNum(dt);
        updateAge(dt);
        //SH_NEW 9_14
        //updatePosition(pos);
        interpolateEnvVars(pos);
        //check for exiting grid
        //SH_NEW: remove this so they don't die when nearshore
        /*if (i3d.isAtGridEdge(pos,tolGridEdge)){
            alive=false;
            active=false;
        }*/
        if (debug) {
            logger.info(toString());
        }
        updateAttributes(); //update the attributes object w/ nmodified values
    }
    
    //WTS_NEW 2012-07-26:{
    //deleted methods calcW(dt) and calcUV(dt)
    
    /**
     * Function to calculate movement rates.
     * 
     * @param dt - time step
     * @return 
     */
    public double[] calcUVW(double[] pos, double dt) {
        //compute vertical velocity
        double w = 0;
        double TL;
        if (fcnVM instanceof wts.models.DisMELS.IBMFunctions.Movement.DielVerticalMigration_FixedDepthRanges) {
            //calculate the vertical movement rate
                    //SH_NEW    
                //w = (Double) fcnVV.calculate(new double[]{dt});
                //Transform length (which is SL) to TL. From T. Hurst
                TL = (length + 0.5169)/0.9315;
                //Calculate swimspeed, ie w (mm/sec.  From T. Hurst
                //w = (0.081221+(0.043168*Math.log10(T)))*Math.pow(TL,1.49652);
                //Make w meters/sec
                //w=w/1000.0;
                //density = w;
                //w = 0.0;
            
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
            double[] ss = DateTimeFunctions.computeSunriseSunset(lon,lat,ModelCalendar.getCalendar().getYearDay());
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
        age        = age+dt/86400;
        ageInStage = ageInStage+dt/86400;
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
        if (fcnMortality instanceof ConstantMortalityRate){
            mortalityRate = (Double)fcnMortality.calculate(null);
        } else 
        if (fcnMortality instanceof InversePowerLawMortalityRate){
          //  mortalityRate = (Double)fcnMortality.calculate(diam);//using egg diameter as covariate for mortality
        }
        double totRate = mortalityRate;
       /* if ((ageInStage>=minStageDuration)) {
            totRate += stageTransRate;
            //apply mortality rate to previous number transitioning and
            //add in new transitioners
            numTrans = numTrans*Math.exp(-dt*mortalityRate/86400)+
                    (stageTransRate/totRate)*number*(1-Math.exp(-dt*totRate/86400));
        }*/
        //if total depth is appropriate for settlement and 
        //indiv is near the bottom, then settle and transform to next stage.
        if ((totalDepth>=minSettlementDepth)&&
                (totalDepth<=maxSettlementDepth)&&
                (depth>(totalDepth+5))) {
            totRate += stageTransRate;
            //apply mortality rate to previous number transitioning and
            //add in new transitioners
            numTrans = numTrans*Math.exp(-dt*mortalityRate/86400)+
                    (stageTransRate/totRate)*number*(1-Math.exp(-dt*totRate/86400));
       
        number = number*Math.exp(-dt*totRate/86400);
        //}: WTS_NEW 2012-07-26
    }
    }
    
    private void updatePosition(double[] pos) {
        totalDepth = i3d.interpolateBathymetricDepth(pos);
        depth      = -i3d.calcZfromK(pos[0],pos[1],pos[2]);
        lat        = i3d.interpolateLat(pos);
        lon        = i3d.interpolateLon(pos);
        gridCellID = ""+Math.round(pos[0])+"_"+Math.round(pos[1]);
        updateTrack();
    }
    
    private void interpolateEnvVars(double[] pos) {
        temperature = i3d.interpolateTemperature(pos);
        salinity    = i3d.interpolateSalinity(pos);
        if (i3d.getPhysicalEnvironment().getField("rho")!=null) rho  = i3d.interpolateValue(pos,"rho");
        else rho = 0.0;
    }

/************************************************************************/
/*	trian(tmin,tmode,tmax)	-SH code 8/2012		*/
/************************************************************************/

/*trian.c - Program returns a random deviate from a triangular       */
/*  distribution.			  	  Oct. 25,1996  -SH  */

/* tmin is the minimum value
   tmode is the mode
   tmax is the maximum value
   tdev is the returned deviate
*/

private double trian(double tmin,double tmode,double tmax) 

{
  int i;
  double tdev=0.0,u,x=0.0;

/*Generate triangular deviate for t(0,1)*/

    u = Math.random();

    if(u<=0.5)  x = Math.sqrt(0.5*u);
    if(u>0.5)   x = 1.0 - Math.sqrt(0.5*(1.0-u));
 
    if(x<0.0)   x = 0.0;
    if(x>1.0)   x = 1.0;

/*Convert to triangular with tmin,tmode,tmax*/

    if(x<=0.5)  tdev = tmin + 2.0*(tmode-tmin)*x;
    if(x>0.5)   tdev = 2.0*tmode-tmax+2.0*(tmax-tmode)*x;

    return tdev;

}

/************************************************************************/
/*				end trian()				*/
/************************************************************************/

    @Override
    public double getStartTime() {
        return startTime;
    }

    @Override
    public void setStartTime(double newTime) {
        startTime = newTime;
        time      = startTime;
        atts.setValue(BenthicJuvStageAttributes.PROP_startTime,startTime);
        atts.setValue(BenthicJuvStageAttributes.PROP_time,time);
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
        atts.setValue(BenthicJuvStageAttributes.PROP_density,density);
        atts.setValue(BenthicJuvStageAttributes.PROP_devStage,devStage);
        atts.setValue(BenthicJuvStageAttributes.PROP_length,length);
        atts.setValue(BenthicJuvStageAttributes.PROP_rho,rho);
        atts.setValue(BenthicJuvStageAttributes.PROP_salinity,salinity);
        atts.setValue(BenthicJuvStageAttributes.PROP_temperature,temperature);    
    }

    /**
     * Updates local variables from the attributes.  
     */
    @Override
    protected void updateVariables() {
        super.updateVariables();
        density     = atts.getValue(BenthicJuvStageAttributes.PROP_density,density);
        devStage    = atts.getValue(BenthicJuvStageAttributes.PROP_devStage,devStage);
        length        = atts.getValue(BenthicJuvStageAttributes.PROP_length,length);
        rho         = atts.getValue(BenthicJuvStageAttributes.PROP_rho,rho);
        salinity    = atts.getValue(BenthicJuvStageAttributes.PROP_salinity,salinity);
        temperature = atts.getValue(BenthicJuvStageAttributes.PROP_temperature,temperature);
     }

}
