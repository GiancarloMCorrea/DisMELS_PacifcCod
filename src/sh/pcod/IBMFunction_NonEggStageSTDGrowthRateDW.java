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
 * IBM function to calculate  size/temperature-dependent growth (STDG) rate in
 * dry weight for non-egg stages using
 *   rate = ((0.454 + 1.610*t - 0.069*t*t)*exp(-6.725*m)+3.705)/100 in g/g/d dry weight
 * where t is temperature in deg C and m is dry weight in micrograms.
 * 
 * From Hurst et al. (2010), equation 4.
 * 
 * @author WilliamStockhausen
 */
@ServiceProviders(value={
    @ServiceProvider(service=IBMGrowthFunctionInterface.class),
    @ServiceProvider(service=IBMFunctionInterface.class)}
)

public class IBMFunction_NonEggStageSTDGrowthRateDW extends AbstractIBMFunction implements IBMGrowthFunctionInterface {
    public static final String DEFAULT_type = "Growth";
    /** user-friendly function name */
    public static final String DEFAULT_name = "Intrinsic growth rate (g/g/d) in dry weight for Pacific cod non-egg stages";
    /** function description */
    public static final String DEFAULT_descr = "Intrinsic growth rate (g/g/d) in dry weight for Pacific cod non-egg stages";
    /** full description */
    public static final String DEFAULT_fullDescr = 
        "\n\t**************************************************************************"+
        "\n\t* This function provides an implementation of the Hurst et al. (2010)"+
        "\n\t* temperature-dependent function for growth in dry weight for Pacific cod Epijuv."+
        "\n\t* "+
        "\n\t* "+
        "\n\t* @author William Stockhausen"+
        "\n\t* "+
        "\n\t* Variables:"+
        "\n\t*      t - Double value of temperature (deg C)"+
        "\n\t*      m - Double value of dry weight (micrograms)"+
        "\n\t* Value:"+
        "\n\t*      r - Double - intrinsic STDG growth rate in dry weight (g/g/d) for egg stages"+
        "\n\t* Calculation:"+
        "\n\t*     r = ((0.454 + 1.610*t - 0.069*t*t)*exp(-6.725*m))/100; (original in %/d)"+
        "\n\t* "+
        "\n\t*  Citation:"+
        "\n\t* Hurst et al. 2010, eq. 4."+
        "\n\t**************************************************************************";
    /** number of settable parameters */
    public static final int numParams = 0;
    /** number of sub-functions */
    public static final int numSubFuncs = 0;
    public IBMFunction_NonEggStageSTDGrowthRateDW(){
        super(numParams,numSubFuncs,DEFAULT_type,DEFAULT_name,DEFAULT_descr,DEFAULT_fullDescr);
    }
    
    @Override
    public Object clone() {
        IBMFunction_NonEggStageSTDGrowthRateDW clone = new IBMFunction_NonEggStageSTDGrowthRateDW();
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
     * Calculates growth rate in dry weight (g/g/d) based on input temperature. 
     * 
     * @param o - Double[] with values 
     *     [0]: in situ temperature in deg C
     *     [1]: dry weight in micrograms
     * 
     * @return Double - growth rate (g/g//d)
     * 
     */
    @Override
    public Object calculate(Object o) {
        Double[] vals = (Double[]) o;
        double t = vals[0];
        double m = vals[1];
        double r = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-6.725*m))/100;// original in %/d
        return (Double) r;
    }
    
}
