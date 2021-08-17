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
 * IBM function to calculate prey-density and temperature-dependent growth (BIOEN) rate in
 * dry weight for non-egg stages.
 * 
 * From Kristiansen et al. (2007).
 * 
 * @author GiancarloMCorrea
 */
@ServiceProviders(value={
    @ServiceProvider(service=IBMGrowthFunctionInterface.class),
    @ServiceProvider(service=IBMFunctionInterface.class)}
)

public class IBMFunction_NonEggStageBIOENGrowthRateDW extends AbstractIBMFunction implements IBMGrowthFunctionInterface {
    public static final String DEFAULT_type = "Growth";
    /** user-friendly function name */
    public static final String DEFAULT_name = "Bioenergetics-Intrinsic growth rate (g/g/d) in dry weight for Pacific cod non-egg stages";
    /** function description */
    public static final String DEFAULT_descr = "Bioenergetics-Intrinsic growth rate (g/g/d) in dry weight for Pacific cod non-egg stages";
    /** full description */
    public static final String DEFAULT_fullDescr = 
        "\n\t**************************************************************************"+
        "\n\t* This function provides an implementation of the Kristiansen et al. (2007)"+
        "\n\t* prey-density/temperature-dependent function for growth in dry weight for"+
        "\n\t* Pcod larvae."+
        "\n\t* "+
        "\n\t* "+
        "\n\t* @author Giancarlo M. Correa"+
        "\n\t* "+
        "\n\t* Variables:"+
        "\n\t*      t - Double value of temperature (deg C)"+
        "\n\t*      m - Double value of dry weight (micrograms)"+
        "\n\t*      OTHERS TO BE DEFINED"+
        "\n\t* Value:"+
        "\n\t*      r - Double - intrinsic BIOEN growth rate in dry weight (g/g/d) for egg stages"+
        "\n\t* Calculation:"+
        "\n\t*     See main paper."+
        "\n\t* "+
        "\n\t*  Citation:"+
        "\n\t* Kristiansen et al. (2007)"+
        "\n\t**************************************************************************";
    /** number of settable parameters */
    public static final int numParams = 0;
    /** number of sub-functions */
    public static final int numSubFuncs = 0;
    public IBMFunction_NonEggStageBIOENGrowthRateDW(){
        super(numParams,numSubFuncs,DEFAULT_type,DEFAULT_name,DEFAULT_descr,DEFAULT_fullDescr);
    }
    
    @Override
    public Object clone() {
        IBMFunction_NonEggStageBIOENGrowthRateDW clone = new IBMFunction_NonEggStageBIOENGrowthRateDW();
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
     *     [2]: dt
     *     [3]: deltaH
     * 
     * @return Double - growth rate (g/g//d)
     * 
     */
    @Override
    public Object calculate(Object o) {
        Double[] vals = (Double[]) o;
        double t = vals[0];
        double m = vals[1];
        double dt = vals[2];
        double deltaH = vals[3];
        double sec2day = 1/86400;
        double meta = dt*2.38E-7*Math.exp(0.088*t)*(Math.pow(m*1000,0.9)*0.001)*deltaH;// mg2ug=1000 here 
        /**
         * Add effect of light
         */
        double assi = 0.8*(1-0.4*Math.exp(-0.002*(m*1000-50)));// mg2ug=1000 here 
        double gr = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-6.725*m)+3.705)/100;// original in %/d
        double g = Math.max(0,(Math.log(gr/100 + 1))*sec2day*dt*deltaH);
        double r = m*(Math.exp(g)-1);
        return (Double) r;
    }
    
}
