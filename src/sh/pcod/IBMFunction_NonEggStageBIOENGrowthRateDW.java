/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package sh.pcod;

import org.openide.util.lookup.ServiceProvider;
import org.openide.util.lookup.ServiceProviders;
import java.util.Arrays;
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
        double attCoeff = vals[10]; // K parameter in Fiksen et al 2002
        // Prey items: should be in this order (by size):
        double eup = vals[11]; // prey item 3
        double ncas = vals[12]; // prey item 3
        double ncao = vals[13]; // prey item 4
        double cop = vals[14]; // prey item 5
        double pCO2 = vals[15]; // pco2 conc
        double ageFromYSL = vals[16];
        double dwmax = vals[17];

        // prey information
        double randomAbunEup = Math.random()*(0.5 + 0.5) - 0.5; // generate random number between -0.5 and 0.5
        double randomAbunNcas = Math.random()*(0.5 + 0.5) - 0.5; // generate random number between -0.5 and 0.5
        double randomAbunNcao = Math.random()*(0.5 + 0.5) - 0.5; // generate random number between -0.5 and 0.5
        double randomAbunCop = Math.random()*(0.5 + 0.5) - 0.5; // generate random number between -0.5 and 0.5

        int npreyitems = 4; // number of prey items
        // include spatial variability:
        double[] zoo_carbon = {eup*Math.exp(randomAbunEup), ncas*Math.exp(randomAbunNcas), ncao*Math.exp(randomAbunNcao), cop*Math.exp(randomAbunCop)}; 
        double[] par_a = {1.38E-8, 2.75E-12, 1E-10, 2.4E-8}; 
        double[] par_b = {2.92, 4.03, 3.56, 2.85}; 
        double[] min_len = {3000, 400, 200, 200}; // minimum length in um
        double[] dlen = {3000, 200, 200, 200}; // size bin in um
        int[] nsizes = {10, 14, 7, 7}; // n size categories 
        // Zooplankton total len vector: (in mm)
        double[] zoolen = {0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30};
        int nallsizes = zoolen.length;
        int[][] zooInd = new int[npreyitems][nallsizes];
        int[] eupInd =  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1}; // for eup
        int[] ncasInd = {0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0}; // for ncas
        int[] ncaoInd = {1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // for ncao
        int[] copInd =  {1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // for cop
        // double[][] prey_len = new double[npreyitems][nallsizes];
        double[][] prey_wgt = new double[npreyitems][nallsizes];
        double[][] prey_area = new double[npreyitems][nallsizes];
        double[][] prey_abun = new double[npreyitems][nallsizes];

        for(int pit=0; pit<npreyitems; pit++) {

            // Zooplankton per len bin:
            double[][] out_zoo = null;
            out_zoo = zooplankton(zoo_carbon[pit], par_a[pit], par_b[pit], min_len[pit], dlen[pit], nsizes[pit]);

            int sind = 0;
            for(int psi=0; psi < nallsizes; psi++) {

                if(pit == 0) zooInd[pit][psi] = eupInd[psi];
                if(pit == 1) zooInd[pit][psi] = ncasInd[psi];
                if(pit == 2) zooInd[pit][psi] = ncaoInd[psi];
                if(pit == 3) zooInd[pit][psi] = copInd[psi];

                if(zooInd[pit][psi] == 1) {
                    // prey_len[pit][psi] = out_zoo[0][sind];
                    prey_wgt[pit][psi] = out_zoo[1][sind];
                    prey_area[pit][psi] = out_zoo[2][sind];
                    prey_abun[pit][psi] = out_zoo[3][sind];
                    sind += 1;
                }

            }

        }

        //calculate pco2 factor:
        double facCO2 = calcCO2(pCO2);

        // EFFECT PCO2: ON METABOLISM
        // Begin function:
        double meta = dtday*2.38e-7*Math.exp(0.088*t)*Math.pow(m,0.9)*(1 + facCO2*0.1); // as in Kristiansen et al 2007. Units: mg/day (without dt). HERE I CHANGED dt FOR dtday 
        double metamax = dtday*2.38e-7*Math.exp(0.088*t)*Math.pow(dwmax,0.9)*(1 + facCO2*0.1); 
        // *(1 + facCO2*0.1)
        // dtday makes more sense 

        if(eb > 0.001) {
            if(std_len > 5.5){
                meta *= 2.5;
                metamax *= 2.5;
            } else {
                meta *= 1.4;
                metamax *= 1.4;
            }
        } 

        double assi = 0.8*(1-0.4*Math.exp(-0.002*(m*1000-50)));// mg2ug=1000 here. No units
        //double r = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-6.725*m))/100;// units: 1/day. This is similar to 'g' (1/day) in TROND. 

        double r = 0;
        r = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-2.225*m))/100;
        double rmax = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-2.225*dwmax))/100;

        // EFFECT PCO2: ON GROWTH
        if(ageFromYSL <= 14) {
            r *= (1 - facCO2*0.1);
            rmax *= (1 - facCO2*0.1);
        }
        if((ageFromYSL > 14) && (ageFromYSL <= 35)) {
            r *= (1 + facCO2*0.15);
            rmax *= (1 + facCO2*0.15);
        }
        // growth1: *(1 - facCO2*0.1)
        // growth2: *(1 + facCO2*0.15)
        double gr_mg = m*(Math.exp(r*dtday) - 1); // Same as TROND
        double gr_mg_max = dwmax*(Math.exp(rmax*dtday) - 1); // Same as TROND

        // START FORAGING PART:
        double contrast = 0.3;
        double em = Math.pow(std_len, 2)/(contrast*0.1*0.2*0.75);
        int dt_num = 100;
        double dt_pca = 0.1;
        double dr = 0.0;
        double pt = 0.0;
        double m_teta = Math.PI/6;
        double var_teta = Math.PI/6;
        double x_star = 0.5*std_len;
        // Values in Fiksen and MacKenzie 2002b:
        double speed_fish = 10.0;
        double speed_prey = 100.0;
        double va = speed_fish*std_len;
        double ke_larvae = 1;
        double beamAttCoeff = attCoeff*3;
        double ke_predator = 1;

        double numing = 0;
        double dening = 0;
        double ing = 0;
        double enc = 0;
        double hand = 0;
        double pca = 0;
        double psa = 0;
        double prey_normal_speed = 0;

        int n_enc = 10;
        double eps = (5.82*1E-9*Math.pow(Math.sqrt(Math.pow(windX,2) + Math.pow(windY,2)), 3))/(depth+0.1); // Equation 1 in MacKenzie and Leggett 1993
        double gape = Math.exp(-3.720 + 1.818*Math.log(std_len) - 0.1219*Math.pow(Math.log(std_len), 2)); // mouth diameter
        double pl_max = 0.08; // max prey len relative to fish len
        double pl_min = 0.03; // min prey len relative to fish len

        Double[] return_vec = new Double[10]; // value to Return should be specified here
        // This is an 2D array of length = 4

        double max_psize = std_len*pl_max; // maximum prey size allowed in diet based on Munk 1997
        double min_psize = std_len*pl_min; // minimum prey size allowed in diet based on Munk 1997
        double sum_numing = 1E-20; // sum ingestion. Very small number to avoid zero error later on.
        double sum_dening = 1E-20; // sum ingestion. Very small number to avoid zero error later on.
        double stomachFullness = 0; // stomach fullness
        double avgRankNum = 0; // to calculate mean rank
        double avgSizeNum = 0; // to calculate mean size

        // Check prey preference based on size:
        double sizePref = 0.055; // prefered ratio size
        double[] ratioLens = new double[nallsizes]; 
        double[] diffLens = new double[nallsizes]; 
        Integer[] rankLens = new Integer[nallsizes]; 
        for(int k=0; k<nallsizes; k++) {
            ratioLens[k] = zoolen[k]/std_len;
            diffLens[k] = Math.abs(ratioLens[k] - sizePref);
        }
        rankLens = rankify(diffLens, nallsizes); // starts at 1, 2, 3, ...
        int i = 0; // loop indicator over lens

        // Begin loop:
        length_loop: for(int itm = 0; itm < nallsizes; itm++) { // 

            if(eb < 1E-15) { // Here is the bug I got. When Eb is extremely small, no run the visual estimation algorithm, then visual = 0
                break length_loop;
            }

            // Find rank:
            int elementToFind = itm + 1; // because it starts at 1
            i = Arrays.asList(rankLens).indexOf(elementToFind);

            if((zoolen[i] > max_psize) || (zoolen[i] < min_psize)) { // only run loop when prey size is smaller than maximum prey size that a larva can capture
                break length_loop;
            } 

            prey_loop: for(int pit=0; pit<npreyitems; pit++) {

                if(prey_abun[pit][i] < 0.00001) { // when very very small prey abundance is present, do not run the prey loop
                    continue prey_loop;
                }

                if(stomachFullness >= 1) {
                    break length_loop;
                } 

                double ier = 0;
                double visual = Math.sqrt(em*contrast*prey_area[pit][i]*(eb/(ke_larvae+eb)));
                double image = prey_area[pit][i];

                double[] getr_out = new double[2];
                getr_out = getr(visual, beamAttCoeff, contrast, image*1E-06, em, ke_larvae, eb, ier); // mm^2 to m^2 = 1000
                //double[] getr_out = {0,10};
                visual = getr_out[1]*1000; // 0 = new ier, 1 = new 'visual' value after getr. m2mm = 1000

                // pca[i] = Math.max(0, Math.min(1, -16.7*(prey_len[i]/std_len) + (3/2)));

                // Capture and approach probabilities:
                double c = 0.5*gape;
                double rs = c + 0.1*std_len;
                double d_crit = 0.264/zoolen[i];
                double w = speed_prey*zoolen[i]; // prey escape velocity
                double capt_pca = 0;
                double capt_psa = 0;
                double travel = 0.43; // s-1. Fiksen and McKenzie 2002
                int max_ats = 3; // max number of ats
                pt = 0;

                if(std_len <= 17) { // Run as Fiksen McKenzie 2002

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

                            // Equation 11 in Fiksen and MacKenzie 2002
                            double capture = (w/va) * (Math.sin(teta)*(rs+c)+(gape/2)*Math.cos(teta));

                            if(capture < (gape*0.5)) {
                                capt_pca = capt_pca + 1;
                            } else {
                                continue enc_loop_part2;
                            }

                        }

                    } // End enc_loop

                    pt = pt/n_enc;

                    psa = Math.min(1, capt_psa/n_enc);
                    pca = Math.min(1, capt_pca/n_enc);

                } // end If std_len <= 17. 

                if(std_len > 17) { //Run as Letcher et al 1996:

                    double par_a_cs = 1.1*std_len/(pl_max*std_len); // pl_max assumed to be 0.09*std_len                                   
                    pca = Math.max(0, 1 - (par_a_cs * (zoolen[i]/std_len)));

                } // end If std_len > 17

                double omega = Math.sqrt(3.615*Math.pow((eps*visual*0.001), 0.667)); //mm2m = 0.001. Equation 11 in MacKenzie Miller 1994. Is it wrong in TROND?
                omega = omega * 1000; // m2mm = 1000. From m/s to mm/s

                // Equation based on Bradley et al 2013, Figure 6:
                prey_normal_speed = zoolen[i]*(1.94*Math.pow(zoolen[i], -1.005)); // this is different from escape velocity
                // Figure 2 in Walton 1992:
                hand = Math.exp(0.264*Math.pow(10, (7.0151*(zoolen[i]/std_len))));
                // See Fiksen and MacKenzie 2002 Equation 1: 
                // EFFECT PCO2: ON PREY ABUNDANCE
                enc = ((0.667*Math.PI*Math.pow(visual,3)*travel + Math.PI*Math.pow(visual,2)*Math.sqrt(Math.pow(prey_normal_speed, 2) + 2*Math.pow(omega,2))*travel*2)*(1*prey_abun[pit][i])*(1 - facCO2*0.1)*1e-6); // tau = 2. MultiplyPrey = 1. ltr2mm3 = 1E-6
                // *(1 - facCO2*0.1)
                // EFFECT PCO2: ON PREY WEIGHT
                // EFFECT PCO2: ON ATTACK SUCCESS
                numing = enc*pca*(1-facCO2*0.1)*prey_wgt[pit][i]*(1-facCO2*0.1)*0.001; // ug2mg = 0.001. 
                // pca: *(1-facCO2*0.1)
                // wgt: *(1-facCO2*0.1)
                dening = enc*hand;

                sum_numing += numing;
                sum_dening += dening;
                // calculate ingestion : Letcher et al 1996
                ing = dt*sum_numing/(1 + sum_dening); 

                avgRankNum += numing*(pit+1); // numing makes more sense than ing
                avgSizeNum += numing*(zoolen[i]);
                // Here calculate suming and check stomach fullness:
                stomachFullness = Math.min(1, (stm_sta + ing/(m*0.06)));

            } // end prey_loop

        } // end of length_loop

        return_vec[0] = gr_mg; // same as TROND
        return_vec[1] = meta; // metabolism
        return_vec[2] = ing;
        return_vec[3] = assi;
        return_vec[4] = stomachFullness;
        return_vec[5] = avgRankNum/sum_numing;
        return_vec[6] = avgSizeNum/sum_numing;
        return_vec[7] = eps;
        return_vec[8] = metamax;
        return_vec[9] = gr_mg_max;
        return return_vec;

    }

    /**
     * Calculates light intensity 
     */
    public static double[] calcLight(double chla, double depth, double bathy) {
        double[] outp = new double[2];
        double attCoef = 0.034 + 0.0518*Math.pow(chla, 0.428) + 0.0363 + 2.833*Math.pow(bathy, -1.079); // Eq A14 in Kearney et al 2020
        double eb_tmp = Math.exp(-1*depth*attCoef);

        // create output:
        outp[0] = attCoef; // K parameter
        outp[1] = eb_tmp; // second part of Eb equation 

        return outp;

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
        sstmp[0] = radfl0; // in W/m^2
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


    public static double[][] zooplankton(double prey_item_totabun_1, double par_a, double par_b, double min_len, double dlen, int nsizes) {

            // size spectrum per prey item:
            // Microzooplankton (MZL): 20-200 um. a=2.63E-6, b=2.23. Stoecker et al 2014, Rodriguez and Mullin 1986
            // Pseudocalanus spp. (Cop): 400 - 1300 um. a=2.4E-8, b=2.85 Liu and Hopcroft 2008
            // Calanus marshallae (NCaS): 900 - 3000 um. a=9.53E-6, b=4.034 Liu and Hopcroft 2007
            // Thysanoessa raschii (EupS): 12000-25000 um. a=7.83E-9, b=3.02 Harding 1977, Becker and Warren 2014
            // Neocalanus spp (NCaO): 400-1400 um. a=1E-10, b = 3.65 Liu and Hopcroft 2006
            // Thysanoessa inermis (EupO): 12000-30000 um. a=1.38E-8, b=2.98. Silva et al 2017, Saunders et al 2013, Becker and Warren 2014
            // Ratio width:length = 0.4

            double[] prey_wgt = new double[nsizes]; // Number of size categories as Daewel et al 2008.
            double[] prey_width = new double[nsizes];
            double[] prey_length = new double[nsizes];  
            double[] sd_pl = new double[nsizes]; 
            double[] sp_pl = new double[nsizes]; 
            double[] m_new = new double[nsizes]; 
            // Return array: 
            double[][] return_array = new double[4][nsizes]; // should be 16?: (number of prey items*4)x(number of len bins)
            // rows 0, 1, 2, 3. (length mm, weight ug, area mm^2, abundance no.ind/L)

            for(int i = 0; i < prey_length.length; i++){ 
                prey_length[i] = min_len + dlen*i; // in um
            }

            // Find weight and width:
            for(int i = 0; i < prey_length.length; i++){ 
                prey_wgt[i] = par_a*Math.pow(prey_length[i], par_b); // in ug
                prey_width[i] = prey_length[i]*0.3/1000; // width = 0.3*length. In mm
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

            for (int i = 0; i < prey_length.length; i++) {
                return_array[0][i] = prey_length[i] / 1000; // THIS IS PREY LENGTH. from um to mm
            }

            for(int i = 0; i < prey_length.length; i++){
                return_array[1][i] = prey_wgt[i]; // THIS IS PREY WEIGHT in ug
            }

            for(int i = 0; i < prey_length.length; i++){
                return_array[2][i]=0.75*(prey_length[i]/1000)*prey_width[i]; // THIS IS PREY AREA: units mm^2
            }

            // it is assumed that 40% is carbon in one individual
            double prey_item_1_ug = (prey_item_totabun_1*1000)*2.5; // COPEPODS: from mg C/m^3 to ug/m^3

            // Split prey items_ug in size categories:
            for(int i = 0; i < prey_length.length; i++){
                return_array[3][i] = sp_pl[i] * (prey_item_1_ug/prey_wgt[i])/1000; // THIS IS ABUNDANCE PER LEN BIN. units: no. ind/L
            }

            return return_array; 

    }

    public static double[] getr(double r, double c, double c0, double ap, double vc, double ke, double eb, double ier) {

          //  r       : start value of r calculated by EASYR
          //  c       : beam attenuation koefficient (m-1)
          //  c0      : prey inherent contrast (d.l.)
          //  ap      : prey area (mm^2)
          //  vc      : parameter characterising visual capacity (d.l.)
          //            this parameter is denoted E' in Aksnes & Utne (1997)
          //  ke      : saturation parameter (uE m-2 s-1). not important units
          //  eb      : background irradiance at depth DEPTH . not important units

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
        double fr2 = Math.log(Math.abs(c0)*ap*vc); // Equation 10 in Aksnes and Giske 1993
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

    // L-W equations: 
    // parameters obtained from get_LW_parameters.R

    public static double getL_fromW(double wgt, double len) {

        double len_out = Math.pow(wgt/1.976E-06, 1/2.974);

        if(len_out < len) len_out = len;

        return len_out;

    }

    // L-W equations: 
    // parameters obtained from get_LW_parameters.R

    public static double getW_fromL(double len) {

        double wgt_out = 1.976E-06*Math.pow(len, 2.974);

        return wgt_out;

    }

    // Mortality function:

    public static double[] TotalMortality(double larval_mm, double eb, double attCoeff, double new_larva_wgt, double stomachFullness, double dwmax) {

        double[] return_mort = new double[5]; // output object

        double larvalShape = 0.2; // Larval width:length ratio
        double contrast = 0.3;
        double em = 5.0E4;
        double visFieldShape = 0.5;
        double fishSwimVel = 0.10;      // Fish cruising velocity (m/s)
        double aPred = 2.77e-6; // same as 0.01 h-1 as Kristiansen et al 2014. Original value in this study = 2.77e-6
        double bPred = -1.3;    // Parameters for purely size-dependent mortality Fiksen et al. 2002. Original value in this study = -1.3
        double starvationMortality = 1e-5; // original value: 7e-6
        double setMort = 1; // 
        double ke_predator = 1;
        double fishDens = 0.0001; // fish density (fish/m^3)
        double deadThreshold = 0.75; // 75% based on Letcher et al 1996
        double m2mm = 1000;
        double beamAttCoeff = attCoeff*3;
        double kval = 5e-6; // for Fiksen et al 2002. Original value: 7e-6

        // Values of fish and larval length are all in meters in this subroutine,
        // which differs from the rest of the routines that are all in mm.
        double larvalWidth = larvalShape*larval_mm;
        double image = larvalWidth*larval_mm;

        double ier = 0;
        double visual = 0.0;

        // All input to getr is either in m (or per m), or in mm (or per mm). Here we use meter (m):
        if(eb < 1E-15) { 
            visual = 0;
        } else {
            double[] getr_out = new double[2];
            getr_out = getr(visual, beamAttCoeff, contrast, image*1E-06, em, ke_predator, eb, ier); // mm^2 to m^2
            visual = getr_out[1]; // 0 = new ier, 1 = new 'visual' value after getr. in m
        }

        // For Fiksen etal 2002:
        double pe = 0.92/(1 + Math.exp((larval_mm-16.3)/4.13)); //Eq 6 in Fiksen et al 2002 

        // Calculate lethal encounter rate with fish setMort is either 0 (off) or 1 (on)
        // double fishMortality = setMort*(visFieldShape*Math.PI*Math.pow(visual,2)*fishSwimVel*fishDens);
        double fishMortality = kval*pe*Math.pow(visual,2); // visual in m
        double invertebrateMortality = setMort*OtherPred(larval_mm, aPred, bPred);
        double starved = AliveOrDead(new_larva_wgt, stomachFullness, deadThreshold, dwmax);
        double mortality = (invertebrateMortality + fishMortality + starved*starvationMortality);

        return_mort[0] = mortality;
        return_mort[1] = starved;
        return_mort[2] = fishMortality*1000000; // to print it in large numbers
        return_mort[3] = invertebrateMortality*1000000; // to print it in large numbers
        return_mort[4] = starved*starvationMortality*1000000; // to print it in large numbers

        return return_mort;

    }


    public static double OtherPred(double larval_mm, double aPred, double bPred){

        // As in Fiksen and Jorgensen 2011:
        double otherPred = aPred*Math.pow(larval_mm, bPred);
        return otherPred;

    }

    public static double AliveOrDead(double new_larva_wgt, double stomachFullness, double deadThreshold, double dwmax){

        // This will NOT work within this code. Do it externally:
        // 1. Calculate maximum growth at a given time using Hurst et al 2010.
        // 2. Compare it with current weight. Then decide if psur = 0.
        double aliveOrDead = 0;

        if ((stomachFullness < 0.01)) {
            aliveOrDead = 1;
        }

        if (dwmax*deadThreshold > new_larva_wgt) {
        //    Give a very high probability of death when belowe 75% (deadThreshold)
            aliveOrDead = 100000;
        }

        return aliveOrDead;

    }

    public static double calcCO2(double co2_val){

        // Resilient:
        // double outVal = Math.min(1, (1/665)*(Math.exp(co2_val*0.005) - 1)); 
        // Non-resilient:
        double outVal = Math.min(1, Math.max(0, -0.5 + (1/1000)*co2_val)); // range: 500 - 1500 uatm (0 to 1)
        // None:
        // double outVal = 0;

        return outVal;

    }

    // Function to print m Maximum elements
    public static Integer[] rankify(double A[], int n)
    {
        // Rank Vector
        Integer[] R = new Integer[n];
     
        // Sweep through all elements in A
        // for each element count the number
        // of less than and equal elements
        // separately in r and s
        for (int i = 0; i < n; i++) {
            int r = 1, s = 1;
             
            for (int j = 0; j < n; j++)
            {
                if (j != i && A[j] < A[i])
                    r += 1;
                     
                if (j != i && A[j] == A[i])
                    s += 1;    
            }
         
        // Use formula to obtain  rank
        R[i] = r + (int)(s - 1) / (int) 2;
     
        }
     
        return R;
         
    }


    // public static double[] calculateJuv(double t, double m, double dt, double dtday, double std_len, double eb, double windX, double windY, double depth, double stm_sta, double attCoeff, double eupo, double eups, double ncas, double ncao, double cop, double mzl) {
    //     // prey information
    //     double npreyitems = 6; // number of prey items
    //     double[] zoo_carbon = {eupo, eups, ncas, ncao, cop, mzl}; 
    //     double[] par_a = {1.38E-8, 7.83E-9, 2.75E-12, 1E-10, 2.4E-8, 2.63E-6}; 
    //     double[] par_b = {2.92, 3.02, 4.03, 3.56, 2.85, 2.23}; 
    //     double[] min_len = {12000, 12000, 900, 400, 400, 20}; // minimum length in um
    //     double[] dlen = {1600, 1100, 180, 90, 80, 16}; // size bin in um

    //     // Begin function:
    //     double meta = dtday*2.38e-7*Math.exp(0.088*t)*Math.pow(m,0.9); // as in Kristiansen et al 2007. Units: mg/day (without dt). HERE I CHANGED dt FOR dtday 
    //     // dtday makes more sense 

    //     if(eb > 0.001) {
    //         if(std_len > 5.5){
    //             meta *= 2.5;
    //         } else {
    //             meta *= 1.4;
    //         }
    //     } 

    //     double assi = 0.8*(1-0.4*Math.exp(-0.002*(m*1000-50)));// mg2ug=1000 here. No units
    //     double r = ((0.454 + 1.610*t - 0.069*t*t)*Math.exp(-6.725*m))/100;// units: 1/day. This is similar to 'g' (1/day) in TROND. 
    //     double gr_mg = m*(Math.exp(r*dtday) - 1); // Same as TROND

    //     // START FORAGING PART:
    //     double contrast = 0.3;
    //     double em = Math.pow(std_len, 2)/(contrast*0.1*0.2*0.75);
    //     double speed_prey = 100.0;
    //     double speed_fish = 10.0;
    //     double va = speed_fish*std_len;
    //     double ke_larvae = 1;
    //     double beamAttCoeff = attCoeff*3;
    //     double ke_predator = 1;

    //     double[] ing = new double[12];
    //     double[] enc = new double[12];
    //     double[] hand = new double[12];
    //     double[] pca = new double[12];
    //     double[] psa = new double[12];
    //     double capt_prob = 0;
    //     double w = 0;

    //     int n_enc = 10;
    //     double eps = (5.82*1e-9*Math.pow(Math.sqrt(Math.pow(windX,2) + Math.pow(windY,2)), 3))/(depth+0.1); // Equation 1 in MacKenzie and Leggett 1993
    //     double gape = Math.exp(-3.720 + 1.818*Math.log(std_len) - 0.1219*Math.pow(Math.log(std_len), 2));

    //     double[] return_vec = new double[5]; // value to Return should be specified here
    //     // This is an 2D array of length = 4

    //     double sum_ing = 0; // sum ingestion 
    //     double stomachFullness = 0; // stomach fullness
    //     double firstCond = 0;
    //     double w_tmp = 0;

    //     prey_item_loop: for(int pit=0; pit<npreyitems; pit++) {

    //         // Zooplankton per len bin:
    //         double[][] out_zoo = null;
    //         out_zoo = zooplankton(zoo_carbon[pit], par_a[pit], par_b[pit], min_len[pit], dlen[pit]);
    //         double[] prey_len = out_zoo[0]; // in mm per len bin
    //         double[] prey_wgt = out_zoo[1]; // in ug per len bin
    //         double[] prey_area = out_zoo[2]; // in mm^2 per len bin
    //         double[] prey_abun = out_zoo[3]; // in no.ind/L per len bin

    //         w_tmp = speed_prey*prey_len[0]; // prey escape velocity
    //         firstCond = Math.max(0, Math.pow(1-((prey_len[0]*w_tmp)/(std_len*va)), 3)); // use capture probability instad of PCA for juveniles

    //         if(firstCond == 0) { // conditional. 

    //             continue prey_item_loop;

    //         } else {

    //         int ii = prey_area.length; // length of prey_area vector.

    //             // Begin loop:
    //             for(int i = ii-1; i>= 0; i--){

    //             // Prey velocity:
    //             w = speed_prey*prey_len[i]; // prey escape velocity

    //                 if(stomachFullness == 1) {

    //                     break prey_item_loop;

    //                 } else {

    //                     double ier = 0;
    //                     double visual = Math.sqrt(em*contrast*prey_area[i]*(eb/(ke_larvae+eb)));
    //                     double image = prey_area[i];

    //                     if(eb < 1E-50) { // Here is the bug I got. When Eb is extremely small, no run the visual estimation algorithm, then visual = 0

    //                         visual = 0;
    //                         ing[i] = 0; // no ingestion in dark environments

    //                     } else {

    //                         double[] getr_out = new double[2];
    //                         getr_out = getr(visual, beamAttCoeff/1000, contrast, image, em, ke_larvae, eb, ier); // m2mm = 1000
    //                         visual = getr_out[1]; // 0 = new ier, 1 = new 'visual' value after getr

    //                         // Figure 2 in Walton 1992. 
    //                         hand[i] = Math.exp(0.264*Math.pow(10, (7.0151*(prey_len[i]/std_len))));
    //                         // See Fall and Fiksen 2019:
    //                         //capt_prob = Math.max(0, Math.pow(1-((prey_len[i]*w)/(std_len*va)), 3)); // use capture probability instad of PCA for juveniles
    //                         capt_prob = Math.max(0, Math.pow(1-((prey_len[i]*w)/(std_len*va)), 3));
    //                         enc[i] = Math.PI*Math.pow(visual*Math.sin(Math.toRadians(30)), 2)*w*(1*prey_abun[i])*1e-6; // units: N prey/s
    //                         // This equation is similar to Ware 1975:
    //                         ing[i] = dt*enc[i]*capt_prob*prey_wgt[i]*0.001/(1 + hand[i]); // ug2mg = 0.001. 

    //                     } // visual else

    //                     sum_ing += ing[i];
    //                     stomachFullness = Math.min(1, (stm_sta + sum_ing/(m*0.06)));

    //                 } // else stomachFullness

    //             } // End prey length loop

    //         } // end 'if' to run loop over lengths based on prey length < gape

    //     } // End prey items loop

    //     return_vec[0] = gr_mg; // same as TROND
    //     return_vec[1] = meta; // metabolism
    //     return_vec[2] = sum_ing;
    //     return_vec[3] = assi;
    //     return_vec[4] = stomachFullness;
    //     return return_vec;

    // }
   
}


    