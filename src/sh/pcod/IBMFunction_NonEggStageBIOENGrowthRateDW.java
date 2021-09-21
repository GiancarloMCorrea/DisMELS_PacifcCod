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
        // double inputs:
        double t = vals[0];
        double m = vals[1];
        double dt = vals[2];
        double dtday = vals[3];
        double std_len = vals[4];
        double eb = vals[5];
        double windX = vals[6];
        double windY = vals[7];
        double depth = vals[8];
        double stm_sta = vals[9]; // stomach 
        double copepod = vals[10]; // prey item 1

        // Begin function:
        double meta = dtday*2.38e-7*Math.exp(0.088*t)*Math.pow(m,0.9); // as in Kristiansen et al 2007. Units: mg/day (without dt). HERE I CHANGED dt FOR dtday 
        // dtday makes more sense 

        if(eb > 0.001) {
            if(std_len > 5.5){
                meta *= 2.5;
            } else {
                meta *= 1.4;
            }
        } 

        double assi = 0.8*(1-0.4*Math.exp(-0.002*(m*1000-50)));// mg2ug=1000 here. No units
        double r = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-6.725*m))/100;// units: 1/day. This is similar to 'g' (1/day) in TROND. 
        double gr_mg = m*(Math.exp(r*dtday) - 1); // Same as TROND

        // START FORAGING PART:
        double contrast = 0.3;
        double em = Math.pow(std_len, 2)/(contrast*0.1*0.2*0.75);
        double[] ing = new double[12];
        double[] enc = new double[12];
        double[] hand = new double[12];
        double[] pca = new double[12];
        double[] psa = new double[12];
        int dt_num = 100;
        double dt_pca = 0.1;
        double dr = 0.0;
        double pt = 0.0;
        double m_teta = Math.PI/6;
        double var_teta = Math.PI/6;
        double x_star = 0.5*std_len;
        double speed_fish = 10.0;
        double speed_prey = 100.0;
        double va = speed_fish*std_len;
        double ke_larvae = 1;
        double attCoeff = 0.18;
        double beamAttCoeff = attCoeff*3;
        double ke_predator = 1;

        int n_enc = 10;
        double eps = (5.82*1e-9*Math.pow(Math.sqrt(Math.pow(windX,2) + Math.pow(windY,2)), 3))/(depth+0.1);
        double gape = Math.exp(-3.720 + 1.818*Math.log(std_len) - 0.1219*Math.pow(Math.log(std_len), 2));

        Double[] return_vec = new Double[4]; // value to Return should be specified here
        // This is an 2D array of length = 4

        // Zooplankton per len bin:
        double[][] out_zoo = null;
        out_zoo = zooplankton(copepod);
        double[] prey_len = out_zoo[0]; // in mm per len bin
        double[] prey_wgt = out_zoo[1]; // in ug per len bin
        double[] prey_area = out_zoo[2]; // in mm^2 per len bin
        double[] prey_abun = out_zoo[3]; // in no.ind/L per len bin

        int ii = prey_area.length; // length of prey_area vector.

        // Begin loop:
        for(int i=0; i<ii; i++){

            double ier = 0;
            double visual = Math.sqrt(em*contrast*prey_area[i]*(eb/(ke_larvae+eb)));
            double image = prey_area[i];

            double[] getr_out = null;
            getr_out = getr(visual, beamAttCoeff/1000, contrast, image, em, ke_larvae, eb, ier); // m2mm = 1000
            visual = getr_out[1]; // 0 = new ier, 1 = new 'visual' value after getr
            pca[i] = Math.max(0, Math.min(1, -16.7*(prey_len[i]/std_len) + (3/2)));

            // Capture and approach probabilities:
            double c = 0.5*gape;
            double rs = c + 0.1*std_len;
            double d_crit = 0.264/prey_len[i];
            double w = speed_prey*prey_len[i];
            double capt_pca = 0;
            double capt_psa = 0;
            double travel = 0.43;
            int max_ats = 3; // max number of ats
            pt = 0;

            // Calculate the probability of approach and capture
            enc_loop: for(int j = 1; j <= n_enc; j++){

                double d = Math.max(visual, c);
                int k_iter_last = 1;
                
                approach_loop: for(int k = 1; k <= dt_num; k++){

                    double v = 0;
                    k_iter_last = k;
                    if(d > rs) {
                        v = (d_crit*2*(Math.pow(d, 4)))/(3*c*((Math.pow(d,2))-(Math.pow(c, 2))));
                        v = dt_pca*Math.min(std_len, v);
                    } else {
                        capt_psa = capt_psa + 1;
                        break approach_loop;
                    }

                    dr = -1*v*(1-(3*c/(2*d))+Math.pow(c,2)/(2*Math.pow(d,3)));
                    d = d + dr;

                }

                pt = pt + k_iter_last*dt_pca;

                // THIS LOOP IS GENERATED BY MYSELF (Giancarlo). IT IS A BETTER WAY TO PROGRAM THIS PART

                enc_loop_part2: for(int j2 = 1; j2 <= max_ats; j2++) {

                    // Generate random number:
                    double r_rand = Math.random();
                    // Apply n_dev function: (begin)
                    double u1 = Math.max(0.00001, r_rand);
                    var_teta = Math.sqrt(-2*Math.log(u1))*Math.cos(2*Math.PI*r_rand);
                    // (end)

                    double teta = (m_teta - var_teta*m_teta);

                    if(teta > Math.PI) { 
                        teta = 2*Math.PI - teta;
                    }

                    teta = Math.abs(teta);

                    if(teta < Math.PI*0.5) {
                        if((gape*0.5/x_star) > Math.tan(teta)) {
                            if((x_star*Math.cos(teta)/w) < ((rs - c + x_star)/va)) {
                                continue enc_loop_part2;
                            } 
                        }
                    }

                    double capture = (w/va) * (Math.sin(teta)*(rs+c)+(gape/2)*Math.cos(teta));

                    if(capture < (gape*0.5)) {
                        capt_pca = capt_pca + 1;
                    } else {
                        continue enc_loop_part2;
                    }

                }

            }

            pt = pt/n_enc;

            psa[i] = Math.min(1, capt_psa/n_enc);
            pca[i] = Math.min(1, capt_pca/n_enc);

            double omega = 1.9*Math.pow((eps*visual*0.001), 0.667); //mm2m = 0.001 
            omega = omega * 1000; // m2mm = 1000. From m/s to mm/s

            hand[i] = 0.264*Math.pow(10, (7.0151*(prey_len[i]/std_len)));
            enc[i] = ((0.667*Math.PI*Math.pow(visual,3)*travel + Math.PI*Math.pow(visual,2)*Math.sqrt(Math.pow(prey_len[i], 2) + 2*Math.pow(omega,2))*travel*2)*(1*prey_abun[i])*1e-6); // tau = 2. MultiplyPrey = 1. ltr2mm3 = 1E-6
            ing[i] = dtday*enc[i]*pca[i]*prey_wgt[i]*0.001/(1 + hand[i]); // ug2mg = 0.001. dt*deltaH (in TROND) = dt (in DisMELS), that is why I deleted deltaH. I used dtday

        }
        
        double sum_ing = 0;
        for (int i = 0; i < ing.length; i++) {  sum_ing += ing[i]; }
        double stomachFullness = Math.min(1, (stm_sta + sum_ing/(m*0.06))); // gut_size= 0.06
        // TODO: check the usage of stomachFullness in TROND. See behavior component

        return_vec[0] = gr_mg; // same as TROND
        return_vec[1] = meta; // metabolism
        return_vec[2] = sum_ing;
        return_vec[3] = assi;
        return return_vec;

    }

    /**
     * Calculates light intensity 
     */
    public static Object calcLight(Object o) {
        Double[] vals = (Double[]) o;
        double chla = vals[0];
        double depth = vals[1];
        double attCoef = 0.1 + 0.054*Math.pow(chla, 2/3) + 0.0088*chla; // Attenuation coefficient. here k0 = 0.1
        double eb_tmp = Math.exp(-1*depth*attCoef);
        return (Double) eb_tmp;
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
        double eepsil = 1e-9;
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
        sstmp[0] = radfl0;
        sstmp[1] = radmax;
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
            slig = ((0.048 - 1.15e-4)/6)*(12+v)+1.15e-4;
        } else if (v >= -18) {
            slig = (((1.15e-4)-1.15e-5)/6)*(18+v)+1.15e-5;
        } else {
            slig = 1.15e-5;
        }

        // create output:
        sstmp2[0] = height;
        sstmp2[1] = slig;

        return sstmp2;

    }


    public static double[][] zooplankton(double prey_item_totabun_1) {

            double minLen_zoo = 140;// in um
            double maxLen_zoo = 1680;// in um
            double[] prey_length = new double[12]; // Number of size categories as Daewel et al 2008. should be 16? 
            double[] sd_pl = new double[12]; // should be 16? 
            double[] sp_pl = new double[12]; // should be 16? 
            double[] m_new = new double[12]; // should be 16? 
            // Return array: 
            double[][] return_array = new double[4][12]; // should be 16?: (number of prey items*4)x(number of len bins)
            // COPEPODS are rows 0, 1, 2, 3. (length mm, weight ug, area mm^2, abundance no.ind/L)

            for(int i = 0; i < prey_length.length; i++){ 
                prey_length[i] = minLen_zoo*(i+1); // in um
            }

            double tm = 0;
            for(int i = 0; i < prey_length.length; i++){
                sd_pl[i] = 695.73 * Math.exp(-0.0083*prey_length[i]); // units: no ind/L
                m_new[i] = Math.exp(2.772 * Math.log(prey_length[i]) - 7.476); // units: ug C
                tm = tm + (m_new[i]*sd_pl[i]);
            }

            // Percentage of total biomass per size category. sum(sp_pl) = 1
            for(int i = 0; i < prey_length.length; i++){
                sp_pl[i] = (sd_pl[i]*m_new[i])/tm;
            }

            // TODO: find reference to change these numbers:
            double[] prey_wgt = {0.1, 0.25, 1.0, 1.51, 3.76, 6.0, 9.0, 15.0, 17.0, 23.13, 30.0, 45.0}; // THIS IS PREY WEIGHT. units: ug. change if n of categories change. 
            double[] prey_width = {0.08, 0.1, 0.1, 0.15, 0.18, 0.2, 0.22, 0.31, 0.35, 0.38, 0.39, 0.40}; // units: mm. change if n of categories change

            for (int i = 0; i < prey_length.length; i++) {
                return_array[0][i] = prey_length[i] / 1000; // THIS IS PREY LENGTH. from um to mm
            }

            for(int i = 0; i < prey_length.length; i++){
                return_array[1][i] = prey_wgt[i]; // THIS IS PREY WEIGHT
            }

            for(int i = 0; i < prey_length.length; i++){
                return_array[2][i]=0.75*prey_length[i]*prey_width[i]; // THIS IS PREY AREA: units mm^2
            }

            // TODO: find reference for this number: it is assumed that 50% is carbon in one individual
            double prey_item_1_ug = (prey_item_totabun_1*1000)*2; // COPEPODS: from mg C/m^3 to ug/m^3

            // Split prey items_ug in size categories:
            for(int i = 0; i < prey_length.length; i++){
                return_array[3][i] = sp_pl[i] * (prey_item_1_ug/prey_wgt[i])/1000; // THIS IS ABUNDANCE PER LEN BIN. units: no. ind/L
            }

            return return_array; // TODO: output other prey items

    }

    public static double[] getr(double r, double c, double c0, double ap, double vc, double ke, double eb, double ier) {

        // output object
        double[] return_r = new double[2]; // output object

        // Run 'easyr' subroutine: (begin)
        double r2 = Math.abs(c0)*ap*vc*(eb/(ke+eb));
        r = Math.sqrt(r2); // new r
        // Run 'easyr' subroutine: (end)

        double eps = 0.0001;
        int iend = 200; // number of iterations
        double tol = r;

        // Run 'deriv' subroutine: (begin)
        double fr2 = Math.log(Math.abs(c0)*ap*vc);
        double fr1 = Math.log(((ke+eb)/eb)*r*r*Math.exp(c*r));
        double f1 = fr1 - fr2;
        double fder = c + 2/r;
        // Run 'deriv' subroutine: (end)

        double tolf = 100*eps;

        for(int i=1; i<=iend; i++){

            if(f1 != 0) {

                if(fder != 0) {
                    double dx = f1/fder;
                    r = r-dx;
                    if(r < 0) {
                        ier = 3; // r out of allowed range (negative)
                        break;
                    }
                    tol = r;
                        // Run 'deriv' subroutine: (begin)
                        fr2 = Math.log(Math.abs(c0)*ap*vc);
                        fr1 = Math.log(((ke+eb)/eb)*r*r*Math.exp(c*r));
                        f1 = fr1 - fr2;
                        fder = c + 2/r;
                        // Run 'deriv' subroutine: (end)
                    tol = eps;
                    double as = Math.abs(r);

                    if((as-1) > 0) {
                        tol = tol*as;
                    }

                    if((Math.abs(dx)-tol) <= 0){
                        if((Math.abs(f1)-tolf) <= 0) {
                            ier = 0; // same as input value. valid r returned
                            break;
                        } else {
                            continue;
                        }
                    } else {
                        continue;
                    }


                } else { // fder == 0
                    ier = 2; // Return in case of zero divisor.
                    break;
                }

            } else { // f1 == 0

                ier = 0; // same as input value. valid r returned
                break;

            }

        }

        // create output:
        return_r[0] = ier;
        return_r[1] = r;

        return return_r;

    }
   
}
