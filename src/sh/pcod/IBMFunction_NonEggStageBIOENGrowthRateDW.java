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
        double std_len = vals[3];
        double meta = dt*2.38E-7*Math.exp(0.088*t)*(Math.pow(m*1000,0.9)*0.001);// mg2ug=1000 here 

        if(std_len > 5.5){
            meta *= 2.5;
        } else {
            meta *= 1.4;
        }

        /**
         * Add effect of light
         */
        double assi = 0.8*(1-0.4*Math.exp(-0.002*(m*1000-50)));// mg2ug=1000 here 
        double r = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-6.725*m))/100;// units: %/day?
        return (Double) r;
    }

    /**
     * Calculates light intensity 
     */
    public static Object calcLight(Object o) {
        Double[] vals = (Double[]) o;
        double chla = vals[0];
        double depth = vals[1];
        double attCoef = 0.1 + 0.054*Math.pow(chla, 2/3) + 0.0088*chla; // Attenuation coefficient. here k0 = 0.1
        double eb = Math.exp(-1*depth*attCoef);
        return (Double) eb;
    }

    public static double[] calcLightQSW(double lat, double time) {
        double[] sstmp = new double[3]; // output object

        double doy = Math.floor(time); // day of year (julian)
        double hour = (time - Math.floor(time))*24; // hour
        //The fractional year, in radians, is
        double gm = ((2*Math.PI)/365)*(doy-1+(hour-12)/24);
        double deg2rad = Math.toRadians(1.0);

        //the solar declination angle (in radians)
        double decl = 0.006918-0.399912*Math.cos(gm)+0.070257*Math.sin(gm)
                              -0.006758*Math.cos(2*gm)+0.000907*Math.sin(2*gm)
                              -0.002697*Math.cos(3*gm)+0.001480*Math.sin(3*gm);

        double sundv=1.00011+0.001280*Math.sin(gm)+0.034221*Math.cos(gm)
                    +0.000077*Math.sin(2*gm)+0.000719*Math.cos(2*gm);

        double sin2=Math.sin(lat*deg2rad)*Math.sin(decl);
        double cos2=Math.cos(lat*deg2rad)*Math.cos(decl);

        // compute 24 hrs mean solar irrradiance at the marine surface layer (unit: w/m^2)
        double pi2 = 8*Math.atan(1);
        double deg = 360/pi2;
        double rad = pi2/360;
        double eepsil = 1E-9;
        double ifrac = 24;
        double fraci = 1/ifrac;
        double absh2o = 0.09; //absorption of water and ozone
        double s0 = 1365; // w/m^2  solar constant
        double cc = 0; // this is clouds = 0 by now. no cloud data.

        double scosz=0;
        double stot=0;
        double radmax=0;

        // create objects for loop:
        double bioday = 0;
        double biohr = 0;
        double hangle = 0;
        double cosz = 0;
        double srad = 0;
        double sdir = 0;
        double sdif = 0;
        double altdeg = 0;
        double cfac = 0;
        double ssurf = 0;

        // begin Loop:
        for(int i=1; i<=ifrac; i++){
            bioday = doy+(i-0.5)*fraci*0.5;
            biohr = bioday*86400;
            biohr = (biohr+43200)%86400;
            hangle = pi2*biohr/86400;
            cosz = Math.max(0,sin2+cos2*Math.cos(hangle));
            scosz = scosz+cosz;
            srad = s0*sundv*cosz;
            sdir = srad*Math.pow(0.7,(Math.min(100,1/(cosz+eepsil))));
            sdif = ((1-absh2o)*srad-sdir)*0.5; 
            altdeg = Math.max(0, Math.asin(Math.min(1, sin2+cos2)))*deg;
            cfac = (1-0.62*cc + 0.0019*altdeg);
            ssurf = (sdir+sdif)*cfac;
            radmax = Math.max(radmax,ssurf);
            stot = stot+ssurf;
        }

        // create outputs:
        double radfl0 = 0;
        double cawdir = 0;

        scosz = scosz*fraci; // 24-hrs mean of  cosz
        radfl0 = stot*fraci; // 24-hrs mean shortw rad in w/m^2
        cawdir = 1-Math.max(0.15, 0.05/(scosz+0.15));

        // create output:
        sstmp[0] = scosz;
        sstmp[1] = radfl0;
        sstmp[2] = cawdir;

        return sstmp;

    }


    public static double[] calcLightSurlig(double lat, double time, double maxLight) {

        double[] sstmp2 = new double[2]; // output object

        double doy = Math.floor(time); // day of year (julian)
        double hour = (time - Math.floor(time))*24; // hour
        double deg2rad = Math.toRadians(1.0);

        double twlight = 5.76;

        double delta = 0.3979*Math.sin( (0.9856*(doy-80) + 1.9171*(Math.sin(0.9856*doy*deg2rad) - 0.98112))*deg2rad );
        double h12 = delta*Math.sin(lat*deg2rad) - Math.pow(1 - Math.pow(delta, 2), 1/2)
                     *(Math.cos(lat*deg2rad))*(Math.cos(15*12*deg2rad));

        double height = delta*Math.sin(lat*deg2rad) - Math.pow(1-Math.pow(delta, 2), 1/2)*Math.cos(lat*deg2rad)*Math.cos(15*hour*deg2rad);
        double v = Math.asin(height*deg2rad);

        double slig = 0;
        // begin conditionals:
        if(v >= 0) {
            slig = maxLight*(height/h12) + twlight;
        } else if (v >= -6) {
            slig = ((twlight - 0.048)/6)*(6+v)+0.048;
        } else if (v >= -12) {
            slig = ((0.048 - 1.15E-4)/6)*(12+v)+1.15E-4;
        } else if (v >= -18) {
            slig = (((1.15E-4)-1.15E-5)/6)*(18+v)+1.15E-5;
        } else {
            slig = 1.15E-5;
        }

        // create output:
        sstmp2[0] = height;
        sstmp2[1] = slig;

        return sstmp2;

    }
   
}
