/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sh.pcod;

import org.openide.util.lookup.ServiceProvider;
import org.openide.util.lookup.ServiceProviders;
import wts.models.DisMELS.framework.IBMFunctions.AbstractIBMFunction;
import wts.models.DisMELS.framework.IBMFunctions.IBMFunctionInterface;
import wts.models.DisMELS.framework.IBMFunctions.IBMGrowthFunctionInterface;

/**
 * IBM function to calculate size/temperature-dependent growth (STDG) rate in 
 * standard length for non-egg stages using
 *   rate = (0.076 + 0.029*t - 0.00002*t*t)/(1-0.059/exp(L^0.0758)) in mm/d
 * where T is temperature (deg C) and L is SL (mm). 
 * 
 * From Hurst et al. (2010) equation 3.
 * 
 * @author William Stockhausen
 */
@ServiceProviders(value={
    @ServiceProvider(service=IBMGrowthFunctionInterface.class),
    @ServiceProvider(service=IBMFunctionInterface.class)}
)

public class IBMFunction_NonEggStageSTDGrowthRateSL extends AbstractIBMFunction implements IBMGrowthFunctionInterface {
    public static final String DEFAULT_type = "Growth";
    /** user-friendly function name */
    public static final String DEFAULT_name = "STDG rate (mm/d) for standard length in Pacific cod non-egg stages";
    /** function description */
    public static final String DEFAULT_descr = "STDG rate (mm/d) for standard length in Pacific cod non-egg stages";
    /** full description */
    public static final String DEFAULT_fullDescr = 
        "\n\t**************************************************************************"+
        "\n\t* This function provides an implementation of the Hurst et al. (2010)"+
        "\n\t* size/temperature-dependent function for growth (STDG) in standard length of non-egg stages."+
        "\n\t* "+
        "\n\t* "+
        "\n\t* @author William Stockhausen"+
        "\n\t* "+
        "\n\t* Variables:"+
        "\n\t*      t - Double value of temperature (deg C)"+
        "\n\t*      L - Double value of standard length (mm)"+
        "\n\t* Value:"+
        "\n\t*      r - Double - growth rate (mm/d)"+
        "\n\t* Calculation:"+
        "\n\t*     r = (0.076 + 0.029*t - 0.00002*t*t)/(1-0.059/exp(L^0.0758))"+
        "\n\t* "+
        "\n\t*  Citation:"+
        "\n\t* Hurst et al. 2010."+
        "\n\t**************************************************************************";
    /** number of settable parameters */
    public static final int numParams = 0;
    /** number of sub-functions */
    public static final int numSubFuncs = 0;
    public IBMFunction_NonEggStageSTDGrowthRateSL(){
        super(numParams,numSubFuncs,DEFAULT_type,DEFAULT_name,DEFAULT_descr,DEFAULT_fullDescr);
    }
    
    @Override
    public Object clone() {
        IBMFunction_NonEggStageSTDGrowthRateSL clone = new IBMFunction_NonEggStageSTDGrowthRateSL();
        clone.setFunctionType(getFunctionType());
        clone.setFunctionName(getFunctionName());
        clone.setDescription(getDescription());
        clone.setFullDescription(getFullDescription());
        return clone;
    }

    @Override
    public boolean setParameterValue(String param,Object value){
        //no parameters to set
        return true;
    }
    
    /**
     * Calculates growth rate in standard length (mm/d) based on input temperature. 
     * 
     * @param o - Double[] with values 
     *      [0]: in situ temperature in deg C.
     *      [1]: current standard length (mm)
     * 
     * @return Double - growth rate (mm/d)
     * 
     */
    @Override
    public Object calculate(Object o) {
        Double[] vals = (Double[]) o;
        double t = vals[0];
        double L = vals[1];
        double r = (0.076 + 0.029*t - 0.00002*t*t)/(1-0.059/Math.exp(Math.pow(L,0.0758)));
        return (Double) r;
    }
    
}
