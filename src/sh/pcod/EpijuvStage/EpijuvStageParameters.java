/*
 * EpijuvStageParameters.java
 *
 * Revisions:
 * 20181011: 1. Removed FCAT_Development and development functions category.
 *           2. Removed EggAscensionRate function as option for vertical movement.
 * 20181015: 1. Removed FCAT_VerticalVelocity and associated functions category;
 *                function for w is hard-wired in calcUVW.
 * 20190722: 1. Added FCAT_HSM IBMFunction category to incorporate habitat suitabiltiy map-type IBMFunctions.
 * 20190723: 1. Added PARAM_minSettlementHSI IBMParameter to define minimum HSI for settlement.
 * 20210205: 1. Added IBMFunction category FCAT_VerticalVelocity back.
 *           2. Added IBMFunction categories FCAT_GrowthSL, FCAT_GrowthDW, FCAT_GrowthTL, FCAT_GrowthWW.
 * 20210208: 1. Added integer flags for IBMFunctions.
 */

package sh.pcod.EpijuvStage;

import java.beans.PropertyChangeSupport;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;
import org.openide.util.lookup.ServiceProvider;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageBIOENGrowthRateDW;
import sh.pcod.IBMFunction_NonEggStageSTDGrowthRateSL;
import wts.models.DisMELS.IBMFunctions.HSMs.HSMFunction_Constant;
import wts.models.DisMELS.IBMFunctions.HSMs.HSMFunction_NetCDF;
import wts.models.DisMELS.IBMFunctions.HSMs.HSMFunction_NetCDF_InMemory;
import wts.models.DisMELS.IBMFunctions.Mortality.ConstantMortalityRate;
import wts.models.DisMELS.IBMFunctions.Mortality.InversePowerLawMortalityRate;
import wts.models.DisMELS.IBMFunctions.Movement.DielVerticalMigration_FixedDepthRanges;
import wts.models.DisMELS.framework.AbstractLHSParameters;
import wts.models.DisMELS.framework.IBMFunctions.IBMFunctionInterface;
import wts.models.DisMELS.framework.IBMFunctions.IBMParameter;
import wts.models.DisMELS.framework.IBMFunctions.IBMParameterBoolean;
import wts.models.DisMELS.framework.IBMFunctions.IBMParameterDouble;
import wts.models.DisMELS.framework.LifeStageParametersInterface;

/**
 * DisMELS class representing parameters for Pacific cod epipelagic juveniles.
 * 
 * This class uses the IBMParameters/IBMFunctions approach to specifying stage-specific parameters.
 * 
 * @author William Stockhausen/Sarah Hinckley
 */
@ServiceProvider(service=LifeStageParametersInterface.class)
public class EpijuvStageParameters extends AbstractLHSParameters {
    
    public static final long serialVersionUID = 1L;
    
    /** the number of IBMParameter objects defined in the class */
    public static final int numParams = 8;
    public static final String PARAM_isSuperIndividual      = "is a super-individual?";
    public static final String PARAM_horizRWP               = "horizontal random walk parameter [m^2]/[s]";
    public static final String PARAM_minStageDuration       = "min stage duration [d]";
    public static final String PARAM_maxStageDuration       = "max stage duration [d]";
    public static final String PARAM_useRandomTransitions   = "use random transitions";
    public static final String PARAM_minSettlementDepth     = "min settlement depth (m)";
    public static final String PARAM_maxSettlementDepth     = "max settlement depth (m)";
    public static final String PARAM_minSettlementHSI       = "min settlement HSI (0-1)";
    
    /** the number of IBMFunction categories defined in the class */
    public static final int numFunctionCats = 8;
    public static final String FCAT_Mortality        = "mortality";
    public static final String FCAT_GrowthSL         = "growth (SL)";
    public static final String FCAT_GrowthDW         = "growth (DW)";
    public static final String FCAT_GrowthTL         = "growth (TL)";
    public static final String FCAT_GrowthWW         = "growth (WW)";
    public static final String FCAT_VerticalMovement = "vertical movement";
    public static final String FCAT_VerticalVelocity = "vertical velocity";    
    public static final String FCAT_HSM              = "habitat suitability";
    
    public static final int FCN_Mortality_ConstantMortalityRate        = 1;
    public static final int FCN_Mortality_InversePowerLawMortalityRate = 2;
    
    public static final int FCN_GrSL_NonEggStageSTDGrowthRate = 1;
    
    public static final int FCN_GrDW_NonEggStageSTDGrowthRate    = 1;
    public static final int FCN_GrDW_NonEggStageBIOENGrowthRate  = 2;
    
    public static final int FCN_VM_DVM_FixedDepthRanges = 1;
    
    public static final int FCN_VV_Epijuv_VerticalSwimmingSpeed = 1;
    
    public static final int FCN_GrTL_Epijuv_GrowthRate = 1;
    
    public static final int FCN_GrWW_Epijuv_GrowthRate = 1;
    
    public static final int FCN_HSM_Constant        = 1;
    public static final int FCN_HSM_NetCDF          = 2;
    public static final int FCN_HSM_NetCDF_InMemory = 3;
    
    private static final Logger logger = Logger.getLogger(EpijuvStageParameters.class.getName());
    
    /** Utility field used by bound properties.  */
    private transient PropertyChangeSupport propertySupport;
    
    /**
     * Creates a new instance of EpijuvStageParameters.
     */
    public EpijuvStageParameters() {
        super("",numParams,numFunctionCats);
        createMapToParameters();
        createMapToPotentialFunctions();
        propertySupport =  new PropertyChangeSupport(this);
    }
    
    /**
     * Creates a new instance of EpijuvStageParameters
     */
    public EpijuvStageParameters(String typeName) {
        super(typeName,numParams,numFunctionCats);
        createMapToParameters();
        createMapToPotentialFunctions();
        propertySupport =  new PropertyChangeSupport(this);
    }
    
    /**
     * This creates the basic parameters mapParams.
     */
    @Override
    protected final void createMapToParameters() {
        String key;
        key = PARAM_isSuperIndividual;    mapParams.put(key,new IBMParameterBoolean(key,key,false));
        key = PARAM_horizRWP;             mapParams.put(key,new IBMParameterDouble(key,key,new Double(0)));
        key = PARAM_minStageDuration;     mapParams.put(key,new IBMParameterDouble(key,key,new Double(0)));
        key = PARAM_maxStageDuration;     mapParams.put(key,new IBMParameterDouble(key,key,new Double(365)));
        key = PARAM_minSettlementDepth;   mapParams.put(key,new IBMParameterDouble(key,key,new Double(365)));
        key = PARAM_maxSettlementDepth;   mapParams.put(key,new IBMParameterDouble(key,key,new Double(365)));
        key = PARAM_minSettlementHSI;     mapParams.put(key,new IBMParameterDouble(key,key,new Double(0)));
        key = PARAM_useRandomTransitions; mapParams.put(key,new IBMParameterBoolean(key,key,false));
    }

    @Override
    protected final void createMapToPotentialFunctions() {
        //create the map from function categories to potential functions in each category
        String cat; 
        Map<String,IBMFunctionInterface> mapOfPotentialFunctions; 
        IBMFunctionInterface ifi;
        
        cat = FCAT_Mortality;  
        mapOfPotentialFunctions = new LinkedHashMap<>(4); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new ConstantMortalityRate(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        ifi = new InversePowerLawMortalityRate(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_GrowthSL; 
        mapOfPotentialFunctions = new LinkedHashMap<>(2); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new IBMFunction_NonEggStageSTDGrowthRateSL(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_GrowthDW; 
        mapOfPotentialFunctions = new LinkedHashMap<>(3); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new IBMFunction_NonEggStageSTDGrowthRateDW(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        ifi = new IBMFunction_NonEggStageBIOENGrowthRateDW(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_GrowthTL; 
        mapOfPotentialFunctions = new LinkedHashMap<>(2); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new IBMFunction_Epijuv_GrowthRateTL();           mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_GrowthWW; 
        mapOfPotentialFunctions = new LinkedHashMap<>(2); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new IBMFunction_Epijuv_GrowthRateWW();           mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_VerticalMovement;  
        mapOfPotentialFunctions = new LinkedHashMap<>(4); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new DielVerticalMigration_FixedDepthRanges(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_VerticalVelocity;  
        mapOfPotentialFunctions = new LinkedHashMap<>(2); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new IBMFunction_Epijuv_VerticalSwimmingSpeed(); mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        
        cat = FCAT_HSM;  
        mapOfPotentialFunctions = new LinkedHashMap<>(4); 
        mapOfPotentialFunctionsByCategory.put(cat,mapOfPotentialFunctions);
        ifi = new HSMFunction_Constant();          mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        ifi = new HSMFunction_NetCDF();            mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
        ifi = new HSMFunction_NetCDF_InMemory();   mapOfPotentialFunctions.put(ifi.getFunctionName(),ifi);
    }
    
    /**
     * Returns a deep copy of the instance.  Values are copied.  
     * Any listeners on 'this' are not(?) copied, so these need to be hooked up.
     * @return - the clone.
     */
    @Override
    public Object clone() {
        EpijuvStageParameters clone = null;
        try {
            clone = (EpijuvStageParameters) super.clone();
            for (String pKey: mapParams.keySet()) {
                clone.setValue(pKey,this.getValue(pKey));
            }
            for (String fcKey: mapOfPotentialFunctionsByCategory.keySet()) {
                Set<String> fKeys = this.getIBMFunctionKeysByCategory(fcKey);
                IBMFunctionInterface sfi = this.getSelectedIBMFunctionForCategory(fcKey);
                for (String fKey: fKeys){
                    IBMFunctionInterface tfi = this.getIBMFunction(fcKey, fKey);
                    IBMFunctionInterface cfi = clone.getIBMFunction(fcKey,fKey);
                    Set<String> pKeys = tfi.getParameterNames();
                    for (String pKey: pKeys) {
                        cfi.setParameterValue(pKey, tfi.getParameter(pKey).getValue());
                    }
                    if (sfi==tfi) clone.setSelectedIBMFunctionForCategory(fcKey, fKey);
                }
            }
            clone.propertySupport = new PropertyChangeSupport(clone);
        } catch (CloneNotSupportedException ex) {
            ex.printStackTrace();
        }
        return clone;
    }

    /**
     * This method is not supported in this implementation.
     *
     * @param strv - array of values (as Strings) used to create the new instance. 
     *              This should be typeName followed by parameter value (as Strings)
     *              in the same order as the keys.
     * @return 
     */
    @Override
    public EpijuvStageParameters createInstance(final String[] strv) {
        throw new UnsupportedOperationException("Not supported.");
    }

    /**
     * Returns a CSV string representation of the parameter values.
     * This method should be overriden by subclasses that add additional parameters, 
     * possibly calling super.getCSV() to get an initial csv string to which 
     * additional field values could be appended.
     * 
     *@return - CSV string parameter values
     */
    @Override
    public String getCSV() {
        String str = typeName;
        for (String key: mapParams.keySet()) str = str+cc+getIBMParameter(key).getValueAsString();
        return str;
    }
                
    /**
     * Returns the comma-delimited string corresponding to the parameters
     * to be used as a header for a csv file.  
     * This should be overriden by subclasses that add additional parameters, 
     * possibly calling super.getCSVHeader() to get an initial header string 
     * to which additional field names could be appended.
     * Use getCSV() to get the string of actual parameter values.
     *
     *@return - String of CSV header names
     */
    @Override
    public String getCSVHeader() {
        String str = "LHS type name";
        for (String key: mapParams.keySet()) str = str+cc+key;
        return str;
    }

    /**
     * Gets the parameter keys.
     * 
     * @return - keys as String array.
     */
    @Override
    public String[] getKeys(){
        String[] strv = new String[mapParams.keySet().size()];
        return mapParams.keySet().toArray(strv);
    }

    /**
     * Sets parameter value identified by the key and fires a property change.
     * @param key   - key identifying attribute to be set
     * @param value - value to set
     */
    @Override
    public void setValue(String key, Object value) {
        if (mapParams.containsKey(key)) {
            IBMParameter p = mapParams.get(key);
            Object old = p.getValue();
            p.setValue(value);
            propertySupport.firePropertyChange(key,old,value);
        }
    }

    /**
     * Adds a PropertyChangeListener to the listener list.
     * @param l The listener to add.
     */
    @Override
    public void addPropertyChangeListener(java.beans.PropertyChangeListener l) {
        propertySupport.addPropertyChangeListener(l);
    }

    /**
     * Removes a PropertyChangeListener from the listener list.
     * @param l The listener to remove.
     */
    @Override
    public void removePropertyChangeListener(java.beans.PropertyChangeListener l) {
        propertySupport.removePropertyChangeListener(l);
    }
}
