/*
 * EggStageAttributes.java
 *
 * Updated 10/11/2018:
 *   Added "attached" as attribute due to changes in DisMELS framework.
 * 2021-02-04: 1. Removed devStage, cop, eup, and nCa as attributes.
 *             2. Changed "diameter" to "SL" and "egg diameter" to "embryo SL".
 * 2021-02-08: 1. Added attributes embryo dry weight and growth rates for SL and DW.
 * 2021-02-11: 1. Added egg stage progression attribute.
 *
 */

package sh.pcod.EggStage;

import java.util.*;
import java.util.logging.Logger;
import org.openide.util.lookup.ServiceProvider;
import wts.models.DisMELS.framework.AbstractLHSAttributes;
import wts.models.DisMELS.framework.IBMAttributes.IBMAttribute;
import wts.models.DisMELS.framework.IBMAttributes.IBMAttributeBoolean;
import wts.models.DisMELS.framework.IBMAttributes.IBMAttributeDouble;

/**
 * DisMELS class representing attributes for Pacific cod eggs.
 * 
 * @author William Stockhausen
 * @author Sarah Hinckley
 */
@ServiceProvider(service=wts.models.DisMELS.framework.LifeStageAttributesInterface.class)
public class EggStageAttributes extends AbstractLHSAttributes   {
    
    /** Number of new attributes defined by this class */
    public static final int numNewAttributes = 10;
    public static final String PROP_attached    = "attached";
    public static final String PROP_stgProg     = "egg stage progression";
    public static final String PROP_SL          = "embryo SL (mm)";
    public static final String PROP_DW          = "embryo DW (mg)";
    //public static final String PROP_stmsta      = "stomach state";
    public static final String PROP_grSL        = "growth rate for SL (mm/d)";
    public static final String PROP_grDW        = "growth rate for DW (1/d)";
    public static final String PROP_density     = "egg density";
    public static final String PROP_temperature = "temperature deg C";
    public static final String PROP_salinity    = "salinity";
    public static final String PROP_rho         = "in situ density";
    
    protected static final Set<String> newKeys = new LinkedHashSet<>((int)(2*numNewAttributes));
    protected static final Set<String> allKeys = new LinkedHashSet<>((int)(2*(numAttributes+numNewAttributes)));
    protected static final Map<String,IBMAttribute> mapAllAttributes = new HashMap<>((int)(2*(numAttributes+numNewAttributes)));
    protected static final String[] aKeys      = new String[numAttributes+numNewAttributes-1];//does not include typeName
    protected static final Class[]  classes    = new Class[numAttributes+numNewAttributes];
    protected static final String[] shortNames = new String[numAttributes+numNewAttributes];
   
    private static final Logger logger = Logger.getLogger(EggStageAttributes.class.getName());
    
    /**
     * This constructor is provided only to facilitate the ServiceProvider functionality.
     * DO NOT USE IT!!
     */
    public EggStageAttributes(){
        super("NULL");
        finishInstantiation();
    }
    
    /**
     * Creates a new attributes instance with type name 'typeName'.
     */
    public EggStageAttributes(String typeName) {
        super(typeName);
        finishInstantiation();
    }
    
    /**
     * Returns a deep copy of the instance.  Values are copied.  
     * Any listeners on 'this' are not(?) copied, so these need to be hooked up.
     * @return - the clone.
     */
    @Override
    public Object clone() {
        EggStageAttributes clone = new EggStageAttributes(typeName);
        for (String key: allKeys) clone.setValue(key,this.getValue(key));
        return clone;
    }

    /**
     * Returns a new instance constructed from the values of the string[].
     * The first value in the string vector must be the type name.
     * Values are set internally by calling setValues(strv) on the new instance.
     * @param strv - vector of values (as Strings) 
     * @return - the new instance
     */
    @Override
    public EggStageAttributes createInstance(final String[] strv) {
        EggStageAttributes atts = new EggStageAttributes(strv[0]);//this sets atts.typeName
        atts.setValues(strv);
        return atts;
    }
    
    private void finishInstantiation(){
        if (newKeys.isEmpty()){
            //set static field information
            mapAllAttributes.putAll(AbstractLHSAttributes.mapAttributes);//add from superclass
            String key;
            key = PROP_attached;   newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeBoolean(key,"attached"));
            key = PROP_stgProg;    newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"stgProg"));
            key = PROP_SL;         newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"SL"));
            key = PROP_DW;         newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"DW"));
            //key = PROP_stmsta;     newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"stmsta"));
            key = PROP_grSL;       newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"grSL"));
            key = PROP_grDW;       newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"grDW"));
            key = PROP_density;    newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"density"));
            key = PROP_temperature;newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"temp"));
            key = PROP_salinity;   newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"sal"));
            key = PROP_rho;        newKeys.add(key); mapAllAttributes.put(key,new IBMAttributeDouble(key,"rho"));
            allKeys.addAll(AbstractLHSAttributes.keys);//add from superclass
            allKeys.addAll(newKeys);//add from this class
            Iterator<String> it = allKeys.iterator();
            int j = 0; it.next();//skip typeName
            while (it.hasNext()) aKeys[j++] = it.next();
        }
        //set instance information
        Map<String,Object> tmpMapValues = new HashMap<>((int)(2*(numNewAttributes+numAttributes)));
        tmpMapValues.putAll(mapValues);//copy from super
        tmpMapValues.put(PROP_attached,   true);
        tmpMapValues.put(PROP_stgProg,    new Double(0));
        tmpMapValues.put(PROP_SL,         new Double(0));
        tmpMapValues.put(PROP_DW,         new Double(0));
        //tmpMapValues.put(PROP_stmsta,     new Double(0));
        tmpMapValues.put(PROP_grSL,       new Double(0));
        tmpMapValues.put(PROP_grDW,       new Double(0));
        tmpMapValues.put(PROP_density,    new Double(0));
        tmpMapValues.put(PROP_temperature,new Double(-1));
        tmpMapValues.put(PROP_salinity,   new Double(-1));
        tmpMapValues.put(PROP_rho,        new Double(-1));
        mapValues = tmpMapValues;//assign to super
    }

    /**
     * Returns the attribute values as an ArrayList (including typeName).
     * 
     * @return 
     */
    @Override
    public ArrayList getArrayList() {
        ArrayList a = super.getArrayList();
        for (String key: newKeys) a.add(getValue(key));
        return a;
    }

    /**
     * Returns the attributes values (not including typeName) as an Object[].
     * 
     * @return 
     */
    @Override
    public Object[] getAttributes() {
        Object[] atts = new Object[numNewAttributes+AbstractLHSAttributes.numAttributes-1];
        int j = 0;
        Iterator<String> it = allKeys.iterator();
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
        String str = super.getCSV();
        Iterator<String> it = newKeys.iterator();
        while (it.hasNext()) str = str+cc+getValueAsString(it.next());
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
        String str = super.getCSVHeader();
        Iterator<String> it = newKeys.iterator();
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
        String str = super.getCSVHeaderShortNames();
        Iterator<String> it = newKeys.iterator();
        while (it.hasNext()) str = str+cc+mapAllAttributes.get(it.next()).shortName;
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
            for (String key: allKeys){
                classes[j++] = mapAllAttributes.get(key).getValueClass();
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
            for (String key: allKeys){
                shortNames[j++] = mapAllAttributes.get(key).shortName;
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
        super.setValues(strv);//set the standard attribute values
        //set the values of the new attributes
        int j = AbstractLHSAttributes.numAttributes;
        try {
            for (String key: newKeys) setValueFromString(key,strv[j++]);
        } catch (java.lang.IndexOutOfBoundsException ex) {
            //@TODO: should throw an exception here that identifies the problem
            String[] aKeys = new String[EggStageAttributes.allKeys.size()];
            aKeys = EggStageAttributes.allKeys.toArray(aKeys);
                String str = "Missing attribute value for "+aKeys[j-1]+".\n"+
                             "Prior values are ";
                for (int i=0;i<(j);i++) str = str+strv[i]+" ";
                javax.swing.JOptionPane.showMessageDialog(
                        null,
                        str,
                        "Error setting attribute values:",
                        javax.swing.JOptionPane.ERROR_MESSAGE);
                throw ex;
        } catch (java.lang.NumberFormatException ex) {
            String[] aKeys = new String[EggStageAttributes.allKeys.size()];
            aKeys = EggStageAttributes.allKeys.toArray(aKeys);
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
        IBMAttribute att = mapAllAttributes.get(key);
        att.setValue(val);
        String str = att.getValueAsString();
        return str;
    }
    
    @Override
    public void setValueFromString(String key, String value) throws NumberFormatException {
        if (!key.equals(PROP_typeName)){
            IBMAttribute att = mapAllAttributes.get(key);
            att.parseValue(value);
            setValue(key,att.getValue());
        }
    }
}
