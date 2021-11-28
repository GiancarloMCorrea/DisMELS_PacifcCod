/**
 * AbstractNonEggStageAttributes.java
 *
 * Updated:
 * 20181011: 1. Added "attached" as attribute due to changes in DisMELS framework
 * 20190716: 1. Added "hsi" as attribute to incorporate habitat suitability index
 *           2. removed PROP_density and PROP_devStage attributes
 * 20190722: 1. Removed "hsi" attribute since most subclasses don't use it
 * 20210205: 1. Added PROP_DW, changed PROP_length to PROP_SL.
 * 20210206: 1. Renamed class to AbstractLarvalAttributes from AbstractNonEggStageAttributes
 *                in conjunction with creation of AbstractJuvenileAttributes.
 * 20210206: 1. Added growth rates in SL and DW as attributes. 
 */

package sh.pcod;

import java.util.*;
import java.util.logging.Logger;
import wts.models.DisMELS.framework.AbstractLHSAttributes;
import wts.models.DisMELS.framework.IBMAttributes.IBMAttribute;
import wts.models.DisMELS.framework.IBMAttributes.IBMAttributeBoolean;
import wts.models.DisMELS.framework.IBMAttributes.IBMAttributeDouble;

/**
 * DisMELS class representing attributes for larval stage Pacific cod classes.
 */
public abstract class AbstractLarvalAttributes extends AbstractLHSAttributes {
    
    /** Number of new attributes defined by this class */
    public static final int numNewAttributes = 24;
    public static final String PROP_attached    = "attached";
    public static final String PROP_SL          = "standard length (mm)";
    public static final String PROP_DW          = "dry weight (mg)";
    public static final String PROP_stmsta      = "stomach state (units)";
    public static final String PROP_psurvival   = "survival probability";
    public static final String PROP_mortfish   = "mortality fish predation";
    public static final String PROP_mortinv   = "mortality invertebrates";
    public static final String PROP_mortstarv   = "mortality starvation";
    public static final String PROP_dwmax   = "max attainable DW";
    public static final String PROP_avgRank   = "average rank in diet";
    public static final String PROP_avgSize   = "average size in diet";
    public static final String PROP_stomachFullness   = "stomach fullness";
    public static final String PROP_pCO2val    = "pCO2 concentration";
    public static final String PROP_grSL        = "growth rate for SL (mm/d)";
    public static final String PROP_grDW        = "growth rate for DW (1/d)";
    public static final String PROP_temperature = "temperature deg C";
    public static final String PROP_salinity    = "salinity";
    public static final String PROP_rho         = "in situ density";
    public static final String PROP_copepod     = "Small copepods mg/m^3 dry wt C";
    public static final String PROP_neocalanus  = "Neocalanoids mg/m^3 dry wt";
    public static final String PROP_euphausiidShelf  = "Euphausiids Shelf mg/m^3 dry wt C";
    public static final String PROP_euphausiid  = "Euphausiids mg/m^3 dry wt C";
    public static final String PROP_neocalanusShelf  = "Neocalanoids Shelf mg/m^3 dry wt";
    public static final String PROP_microzoo  = "Microzooplankton mg/m^3 dry wt C";

    /** these fields HIDE static fields from superclass and should incorporate ALL information from superclasses */
    protected static final int numAttributes = AbstractLHSAttributes.numAttributes+numNewAttributes;
    protected static final Set<String> keys = new LinkedHashSet<>(2*numAttributes);
    protected static final Map<String,IBMAttribute> mapAttributes = new HashMap<>(2*numAttributes);
    protected static final String[] aKeys      = new String[numAttributes-1];//does not include typeName
    protected static final Class[]  classes    = new Class[numAttributes];
    protected static final String[] shortNames = new String[numAttributes];
   
    private static final Logger logger = Logger.getLogger(AbstractLarvalAttributes.class.getName());
    
    /**
     * This constructor is provided only to facilitate the ServiceProvider functionality.
     * DO NOT USE IT!!
     */
    protected AbstractLarvalAttributes(){
        super("NULL");
        finishInstantiation();
    }
    
    /**
     * Creates a new attributes instance with type name 'typeName'.
     */
    protected AbstractLarvalAttributes(String typeName) {
        super(typeName);
        finishInstantiation();
    }
    
    /**
     * This method adds default values for the new attributes to the superclass field "mapValues".
     * 
     * When the first instance of this class is created, this method also fills in 
     * the static fields "keys" and "mapAttributes" with keys and attributes from the 
     * superclass and from this class.
     */
    private void finishInstantiation(){
        if (keys.isEmpty()){
            //set static field information
            keys.addAll(AbstractLHSAttributes.keys);//add from superclass
            mapAttributes.putAll(AbstractLHSAttributes.mapAttributes);//add from superclass
            String key;
            key = PROP_attached;   keys.add(key); mapAttributes.put(key,new IBMAttributeBoolean(key,"attached"));
            key = PROP_SL;         keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"SL"));
            key = PROP_DW;         keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"DW"));
            key = PROP_stmsta;     keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"stmsta"));
            key = PROP_psurvival;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"psurvival"));
            key = PROP_mortfish;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"mortfish"));
            key = PROP_mortinv;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"mortinv"));
            key = PROP_mortstarv;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"mortstarv"));
            key = PROP_dwmax;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"dwmax"));
            key = PROP_avgRank;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"avgRank"));
            key = PROP_avgSize;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"avgSize")); 
            key = PROP_stomachFullness;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"stomachFullness"));                   
            key = PROP_pCO2val;  keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"pCO2val"));
            key = PROP_grSL;       keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"grSL"));
            key = PROP_grDW;       keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"grDW"));
            key = PROP_temperature;keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"temp"));
            key = PROP_salinity;   keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"sal"));
            key = PROP_rho;        keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"rho"));
            key = PROP_copepod;    keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"copepod"));
            key = PROP_euphausiid; keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"euphausiid"));
            key = PROP_euphausiidShelf; keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"euphausiidShelf"));
            key = PROP_neocalanus; keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"neocalanus"));
            key = PROP_neocalanusShelf; keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"neocalanusShelf"));
            key = PROP_microzoo; keys.add(key); mapAttributes.put(key,new IBMAttributeDouble(key,"microzoo"));

            Iterator<String> it = keys.iterator();
            int j = 0; it.next();//skip typeName
            while (it.hasNext()) aKeys[j++] = it.next();
        }
        //set instance information
        Map<String,Object> tmpMapValues = new HashMap<>(2*numAttributes);
        tmpMapValues.putAll(mapValues);//copy from super
        tmpMapValues.put(PROP_attached,   false);
        tmpMapValues.put(PROP_SL,         new Double(0));
        tmpMapValues.put(PROP_DW,         new Double(0));
        tmpMapValues.put(PROP_stmsta,     new Double(0));
        tmpMapValues.put(PROP_psurvival,  new Double(1));
        tmpMapValues.put(PROP_mortfish,  new Double(1));
        tmpMapValues.put(PROP_mortinv,  new Double(1));
        tmpMapValues.put(PROP_mortstarv,  new Double(1));
        tmpMapValues.put(PROP_dwmax,  new Double(1));
        tmpMapValues.put(PROP_avgRank,  new Double(1));
        tmpMapValues.put(PROP_avgSize,  new Double(1));
        tmpMapValues.put(PROP_stomachFullness,  new Double(1));
        tmpMapValues.put(PROP_pCO2val,    new Double(0));
        tmpMapValues.put(PROP_grSL,       new Double(0));
        tmpMapValues.put(PROP_grDW,       new Double(0));
        tmpMapValues.put(PROP_temperature,new Double(-1));
        tmpMapValues.put(PROP_salinity,   new Double(-1));
        tmpMapValues.put(PROP_rho,        new Double(-1));
        tmpMapValues.put(PROP_copepod,    new Double(-1));
        tmpMapValues.put(PROP_euphausiid, new Double(-1));
        tmpMapValues.put(PROP_euphausiidShelf, new Double(-1));
        tmpMapValues.put(PROP_neocalanus, new Double(-1));
        tmpMapValues.put(PROP_neocalanusShelf, new Double(-1));
        tmpMapValues.put(PROP_microzoo, new Double(-1));
        mapValues = tmpMapValues;//assign to super
    }

    /**
     * Returns the attribute values as an ArrayList (including typeName).
     * 
     * @return 
     */
    @Override
    public ArrayList getArrayList() {
        ArrayList a = new ArrayList(keys.size());
        a.add(typeName);
        Iterator<String> it = keys.iterator();
        it.next();//skip PROP_typeName
        while (it.hasNext()) a.add(getValue(it.next()));
        return a;
    }

    /**
     * Returns the attributes values (not including typeName) as an Object[].
     * 
     * @return 
     */
    @Override
    public Object[] getAttributes() {
        Object[] atts = new Object[numAttributes-1];
        int j = 0;
        Iterator<String> it = keys.iterator();
        it.next();//skip PROP_typeName
        while (it.hasNext()) atts[j++] = getValue(it.next()); 
        return atts;
    }
    
   /**
     * Returns a CSV string representation of the attribute values.
     * 
     *@return - CSV string attribute values
     */
    @Override
    public String getCSV() {
        String str = typeName;
        Iterator<String> it = keys.iterator();
        it.next();//skip typeName
        while (it.hasNext()) {
            String key = it.next();
            str = str+cc+getValueAsString(key);
        }
        return str;
    }
                
    /**
     * Returns the comma-delimited string corresponding to the attributes
     * to be used as a header for a csv file.  
     * Use getCSV() to get the string of actual attribute values.
     *
     *@return - String of CSV header names
     */
    @Override
    public String getCSVHeader() {
        Iterator<String> it = keys.iterator();
        String str = it.next();//typeName
        while (it.hasNext()) str = str+cc+it.next();
        return str;
    }
                
    /**
     * Returns the comma-delimited string corresponding to the attributes
     * to be used as a header for a csv file.  
     *
     *@return - String of CSV header names (short style)
     */
    @Override
    public String getCSVHeaderShortNames() {
        Iterator<String> it = keys.iterator();
        String str = mapAttributes.get(it.next()).shortName;//this is "typeName"
        while (it.hasNext())  str = str+cc+mapAttributes.get(it.next()).shortName;
        return str;
    }
    
    /**
     * Returns Class types for all attributes (including typeName) as a Class[]
     * in the order the allKeys are defined.
     * 
     * @return 
     */
    @Override
    public Class[] getClasses() {
        if (classes[0]==null){
            int j = 0;
            for (String key: keys){
                classes[j++] = mapAttributes.get(key).getValueClass();
            }
        }
        return classes;
    }

    /**
     * Returns keys for all attributes excluding typeName as a String[]
     * in the order the keys are defined.
     * 
     * @return 
     */
    @Override
    public String[] getKeys() {        
        return aKeys;
    }

    /**
     * Returns short names for all attributes (including typeName) as a String[]
     * in the order the allKeys are defined.
     * 
     * @return 
     */
    @Override
    public String[] getShortNames() {
        if (shortNames[0]==null){
            int j = 0;
            for (String key: keys){
                shortNames[j++] = mapAttributes.get(key).shortName;
            }
        }
        return shortNames;
    }
    
    /**
     * Sets attribute values to those of input String[].
     * @param strv - String[] of attribute values.
     */
    @Override
    public void setValues(final String[] strv) {
        int j = 1;
        try {
            Iterator<String> it = keys.iterator();
            it.next();//skip typeName
            while (it.hasNext()) setValueFromString(it.next(),strv[j++]);
        } catch (java.lang.IndexOutOfBoundsException ex) {
            //@TODO: should throw an exception here that identifies the problem
            String[] aKeys = new String[keys.size()];
            aKeys = keys.toArray(aKeys);
                String str = "Missing attribute value for "+aKeys[j]+".\n"+
                             "Prior values are ";
                for (int i=0;i<(j);i++) str = str+strv[i]+" ";
                javax.swing.JOptionPane.showMessageDialog(
                        null,
                        str,
                        "Error setting attribute values:",
                        javax.swing.JOptionPane.ERROR_MESSAGE);
                throw ex;
        } catch (java.lang.NumberFormatException ex) {
            String[] aKeys = new String[keys.size()];
            aKeys = keys.toArray(aKeys);
            String str = "Bad attribute value for "+aKeys[j-2]+".\n"+
                         "Value was '"+strv[j-1]+"'.\n"+
                         "Entry was '";
            try {
                for (int i=0;i<(strv.length-1);i++) {
                    if ((strv[i]!=null)&&(!strv[i].isEmpty())) {
                        str = str+strv[i]+", ";
                    } else {
                        str = str+"<missing_value>, ";
                    }
                }
                if ((strv[strv.length-1]!=null)&&(!strv[strv.length-1].isEmpty())) {
                    str = str+strv[strv.length-1]+"'.";
                } else {
                    str = str+"<missing_value>'.";
                }
            }  catch (java.lang.IndexOutOfBoundsException ex1) {
                //do nothing
            }
            javax.swing.JOptionPane.showMessageDialog(
                    null,
                    str,
                    "Error setting attribute values:",
                    javax.swing.JOptionPane.ERROR_MESSAGE);
            throw ex;
        }
    }
    
    @Override
    public String getValueAsString(String key){
        Object val = getValue(key);
        IBMAttribute att = mapAttributes.get(key);
        att.setValue(val);
        String str = att.getValueAsString();
        return str;
    }
    
    @Override
    public void setValueFromString(String key, String value) throws NumberFormatException {
        if (!key.equals(PROP_typeName)){
            IBMAttribute att = mapAttributes.get(key);
            att.parseValue(value);
            setValue(key,att.getValue());
        }
    }
}
