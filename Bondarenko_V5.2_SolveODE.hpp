#include <stdio.h>
#include <math.h>
#include "MCutil.hpp"

/*#include "../../sundials-2.3.0/include/sundials/sundials_types.h"
#include "../../sundials-2.3.0/include/cvode/cvode.h"
#include "../../sundials-2.3.0/include/cvode/cvode_dense.h"
#include "../../sundials-2.3.0/include/cvode/cvode_diag.h"
#include "../../sundials-2.3.0/include/nvector/nvector_serial.h"
#include "../../sundials-2.3.0/include/sundials/sundials_dense.h"

*/

#include "../../sundials/include/sundials/sundials_types.h"
#include "../../sundials/include/cvode/cvode.h"
#include "../../sundials/include/cvode/cvode_dense.h"
#include "../../sundials/include/cvode/cvode_diag.h"
#include "../../sundials/include/nvector/nvector_serial.h"
#include "../../sundials/include/sundials/sundials_dense.h"











#define ABSTOL 1.0E-06
#define RELTOL 1.0E-04
static int check_flag(void *flagvalue, char *funcname, int opt);
static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__);
class Solveode
{
	//PARAMETERS
private:
        int protocolo; //define qual protocolo deve ser utilizado
	int it_countx;
        double v_step;     //voltage step
	double time; 	 // millisecond
	double i_stim; 	 // picoA_per_picoF
	double Cm; 	 // microF_per_cm2
	double Acap; 	 // cm2
	double Vmyo; 	 // microlitre
	double F; 	 // coulomb_per_millimole
	double VJSR; 	 // microlitre
	double Vss; 	 // microlitre
	double VNSR; 	 // microlitre
	double CMDN_tot; 	 // micromolar
	double Km_CMDN; 	 // micromolar
	double CSQN_tot; 	 // micromolar
	double Km_CSQN; 	 // micromolar
	double v1; 	 // per_millisecond
	double tau_tr; 	 // millisecond
	double tau_xfer; 	 // millisecond
	double v2; 	 // per_millisecond
	double v3; 	 // micromolar_per_millisecond
	double Km_up; 	 // micromolar
	double k_plus_htrpn; 	 // per_micromolar_millisecond
	double HTRPN_tot; 	 // micromolar
	double k_plus_ltrpn; 	 // per_micromolar_millisecond
	double LTRPN_tot; 	 // micromolar
	double k_minus_htrpn; 	 // per_millisecond
	double k_minus_ltrpn; 	 // per_millisecond
	double i_CaL_max; 	 // picoA_per_picoF
	double k_plus_a; 	 // micromolar4_per_millisecond
	double n; 	 // dimensionless
	double k_minus_b; 	 // per_millisecond
	double k_minus_c; 	 // per_millisecond
	double k_minus_a; 	 // per_millisecond
	double k_plus_b; 	 // micromolar3_per_millisecond
	double m; 	 // dimensionless
	double k_plus_c; 	 // per_millisecond
	double g_CaL; 	 // milliS_per_microF
	double E_CaL; 	 // millivolt
	double Kpcb; 	 // per_millisecond
	double Kpc_max; 	 // per_millisecond
	double Kpc_half; 	 // micromolar
	double i_pCa_max; 	 // picoA_per_picoF
	double Km_pCa; 	 // micromolar
	double k_NaCa; 	 // picoA_per_picoF
	double K_mNa; 	 // micromolar
	double Nao; 	 // micromolar
	double K_mCa; 	 // micromolar
	double Cao; 	 // micromolar
	double k_sat; 	 // dimensionless
	double eta; 	 // dimensionless
	double R; 	 // joule_per_mole_kelvin
	double T; 	 // kelvin
	double g_Cab; 	 // milliS_per_microF
	double g_Na; 	 // milliS_per_microF
	double Ko; 	 // micromolar
	double g_Nab; 	 // milliS_per_microF
	double g_Kto_f; 	 // milliS_per_microF
	double g_Kto_s; 	 // milliS_per_microF
	double g_Ks; 	 // milliS_per_microF
	double g_Kur; 	 // milliS_per_microF
	double g_Kss; 	 // milliS_per_microF
	double g_Kr; 	 // milliS_per_microF
	double kf; 	 // per_millisecond
	double kb; 	 // per_millisecond
	double i_NaK_max; 	 // picoA_per_picoF
	double Km_Nai; 	 // micromolar
	double Km_Ko; 	 // micromolar
	double g_ClCa; 	 // milliS_per_microF
	double Km_Cl; 	 // micromolar
	double E_Cl; 	 // millivolt
	double IstimAmplitude; 	 // picoA_per_picoF
	double IstimStart; 	 // millisecond
	double IstimEnd; 	 // millisecond
	double IstimPeriod; 	 // millisecond
	double IstimPulseDuration; 	 // millisecond
	double dtime, *time_vec__;
	double time_new;

	//DEPENDENT VARIABLES
	double *V;
	double V_new_, V_old_, V_ini_;
	double *Cai;
	double Cai_new_, Cai_old_, Cai_ini_;
	double *Cass;
	double Cass_new_, Cass_old_, Cass_ini_;
	double *CaJSR;
	double CaJSR_new_, CaJSR_old_, CaJSR_ini_;
	double *CaNSR;
	double CaNSR_new_, CaNSR_old_, CaNSR_ini_;
	double *P_RyR;
	double P_RyR_new_, P_RyR_old_, P_RyR_ini_;
	double *LTRPN_Ca;
	double LTRPN_Ca_new_, LTRPN_Ca_old_, LTRPN_Ca_ini_;
	double *HTRPN_Ca;
	double HTRPN_Ca_new_, HTRPN_Ca_old_, HTRPN_Ca_ini_;
	double *P_O1;
	double P_O1_new_, P_O1_old_, P_O1_ini_;
	double *P_O2;
	double P_O2_new_, P_O2_old_, P_O2_ini_;
	double *P_C2;
	double P_C2_new_, P_C2_old_, P_C2_ini_;
	double *O;
	double O_new_, O_old_, O_ini_;
	double *C2;
	double C2_new_, C2_old_, C2_ini_;
	double *C3;
	double C3_new_, C3_old_, C3_ini_;
	double *C4;
	double C4_new_, C4_old_, C4_ini_;
	double *I1;
	double I1_new_, I1_old_, I1_ini_;
	double *I2;
	double I2_new_, I2_old_, I2_ini_;
	double *I3;
	double I3_new_, I3_old_, I3_ini_;
	double *Nai;
	double Nai_new_, Nai_old_, Nai_ini_;
	double *C_Na2;
	double C_Na2_new_, C_Na2_old_, C_Na2_ini_;
	double *C_Na1;
	double C_Na1_new_, C_Na1_old_, C_Na1_ini_;
	double *O_Na;
	double O_Na_new_, O_Na_old_, O_Na_ini_;
	double *IF_Na;
	double IF_Na_new_, IF_Na_old_, IF_Na_ini_;
	double *I1_Na;
	double I1_Na_new_, I1_Na_old_, I1_Na_ini_;
	double *I2_Na;
	double I2_Na_new_, I2_Na_old_, I2_Na_ini_;
	double *IC_Na2;
	double IC_Na2_new_, IC_Na2_old_, IC_Na2_ini_;
	double *IC_Na3;
	double IC_Na3_new_, IC_Na3_old_, IC_Na3_ini_;
	double *Ki;
	double Ki_new_, Ki_old_, Ki_ini_;
	double *ato_f;
	double ato_f_new_, ato_f_old_, ato_f_ini_;
	double *ito_f;
	double ito_f_new_, ito_f_old_, ito_f_ini_;
	double *ato_s;
	double ato_s_new_, ato_s_old_, ato_s_ini_;
	double *ito_s;
	double ito_s_new_, ito_s_old_, ito_s_ini_;
	double *nKs;
	double nKs_new_, nKs_old_, nKs_ini_;
	double *aur;
	double aur_new_, aur_old_, aur_ini_;
	double *iur;
	double iur_new_, iur_old_, iur_ini_;
	double *aKss;
	double aKss_new_, aKss_old_, aKss_ini_;
	double *iKss;
	double iKss_new_, iKss_old_, iKss_ini_;
	double *C_K2;
	double C_K2_new_, C_K2_old_, C_K2_ini_;
	double *C_K1;
	double C_K1_new_, C_K1_old_, C_K1_ini_;
	double *O_K;
	double O_K_new_, O_K_old_, O_K_ini_;
	double *I_K;
	double I_K_new_, I_K_old_, I_K_ini_;
        double *istim;
         
	// CVODE VARIABLES
	realtype reltol__;
	void *cvode_mem_cvode__;
	N_Vector dependent_variable__, abstol__;
	int flag__, flagr__;
	double *depvar__;

        //currents

        double *I_Ktof;
        double *I_Kur;
        double *I_Kss;
        double *I_Ks;
        double *I_Kr;
        double *I_K1;
        double *I_CaL;
        double *I_Cab;
        double *I_Na;
        double *I_NaK;
        double *I_NaCa;
        double *I_Nab;
        double *I_pCa;
        double *I_ClCa;
        double *I_stim;

public:
	Solveode(double val_abstol__ = ABSTOL, double val_reltol__= RELTOL, int method__ = CV_BDF);
	~Solveode();
	int setVariables(int, double);
	int setParameters(int, double);
	int setFreeVariable(double);
	void setAbstol(int index, double value);
	void setReltol(double value);
	void reInitCVODE();
	double getVariables(int);
	double getParameters(int);
	double getFreeVariable();
	Variables get_Parameters();
	Variables get_Variables();
	Variables get_FreeVariable();
	void setParametersFromFile(char*);
	void setVariablesFromFile(char*);
	void setFreeVariableFromFile(char*);
	double* solveDiff();
	void solveCVODE(int firstcall__ = 0, int steps__ = 0);
	double solve(int firstcall__ = 0, int num_iterations = 0, int num_results__ = 0);
	double* getSolution(int indVariable);
	double* getIndependentVar();
	double solveToFile(char *filename, char *fileaccess = "", int firstcall__ = 0, int num_iterations__ = 0, int num_results__ = 0);
	void solveCVODEToFile(char *filename, char *fileaccess = "", int firstcall__ = 0, int steps__ = 0);
private:
	inline double calc_Bi();
	inline double calc_Bss();
	inline double calc_BJSR();
	inline double calc_J_rel();
	inline double calc_J_tr();
	inline double calc_J_xfer();
	inline double calc_J_leak();
	inline double calc_J_up();
	inline double calc_J_trpn();
	inline double calc_P_C1();
	inline double calc_i_CaL();
	inline double calc_C1();
	inline double calc_alpha();
	inline double calc_beta();
	inline double calc_gamma();
	inline double calc_Kpcf();
	inline double calc_i_pCa();
	inline double calc_i_NaCa();
	inline double calc_i_Cab();
	inline double calc_E_CaN();
	inline double calc_i_Na();
	inline double calc_E_Na();
	inline double calc_C_Na3();
	inline double calc_alpha_Na11();
	inline double calc_alpha_Na12();
	inline double calc_alpha_Na13();
	inline double calc_beta_Na11();
	inline double calc_beta_Na12();
	inline double calc_beta_Na13();
	inline double calc_alpha_Na3();
	inline double calc_beta_Na3();
	inline double calc_alpha_Na2();
	inline double calc_beta_Na2();
	inline double calc_alpha_Na4();
	inline double calc_beta_Na4();
	inline double calc_alpha_Na5();
	inline double calc_beta_Na5();
	inline double calc_i_Nab();
	inline double calc_i_Kto_f();
	inline double calc_E_K();
	inline double calc_alpha_a();
	inline double calc_beta_a();
	inline double calc_alpha_i();
	inline double calc_beta_i();
	inline double calc_i_Kto_s();
	inline double calc_ass();
	inline double calc_iss();
	inline double calc_tau_ta_s();
	inline double calc_tau_ti_s();
	inline double calc_i_K1();
	inline double calc_i_Ks();
	inline double calc_alpha_n();
	inline double calc_beta_n();
	inline double calc_i_Kur();
	inline double calc_tau_aur();
	inline double calc_tau_iur();
	inline double calc_i_Kss();
	inline double calc_tau_Kss();
	inline double calc_i_Kr();
	inline double calc_C_K0();
	inline double calc_alpha_a0();
	inline double calc_beta_a0();
	inline double calc_alpha_a1();
	inline double calc_beta_a1();
	inline double calc_i_NaK();
	inline double calc_f_NaK();
	inline double calc_sigma();
	inline double calc_i_ClCa();
	inline double calc_O_ClCa();
	inline double calc_Istim();
	inline double calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current();
	inline double calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current();
	inline double ifnumber_0();
        inline double ifnumber_1();
        inline double ifnumber_2();
};

#define AGOS_NAN (0.0/0.0)
#define AGOS_INF (1.0/0.0)
#define __agos_xor(a,b) (!(a && b) && (a || b))
float __agos_factorial(int);

	//CONSTRUCTOR
	Solveode::Solveode(double val_abstol__, double val_reltol__, int method__)
	{
		time = 0.000000e+00;
		i_stim = 0.000000e+00;
		Cm = 1.000000e+00;
		Acap = 1.534000e-04;
		Vmyo = 2.584000e-05;
		F = 9.650000e+01;
		VJSR = 1.200000e-07;
		Vss = 1.485000e-09;
		VNSR = 2.098000e-06;
		CMDN_tot = 5.000000e+01;
		Km_CMDN = 2.380000e-01;
		CSQN_tot = 1.500000e+04;
		Km_CSQN = 8.000000e+02;
		v1 = 4.500000e+00;
		tau_tr = 2.000000e+01;
		tau_xfer = 8.000000e+00;
		v2 = 1.740000e-05;
		v3 = 4.500000e-01;
		Km_up = 5.000000e-01;
		k_plus_htrpn = 2.370000e-03;
		HTRPN_tot = 1.400000e+02;
		k_plus_ltrpn = 3.270000e-02;
		LTRPN_tot = 7.000000e+01;
		k_minus_htrpn = 3.200000e-05;
		k_minus_ltrpn = 1.960000e-02;
		i_CaL_max = 7.000000e+00;
		k_plus_a = 6.075000e-03;
		n = 4.000000e+00;
		k_minus_b = 9.650000e-01;
		k_minus_c = 8.000000e-04;
		k_minus_a = 7.125000e-02;
		k_plus_b = 4.050000e-03;
		m = 3.000000e+00;
		k_plus_c = 9.000000e-03;
//		g_CaL = 1.729000e-01;
		E_CaL = 6.300000e+01;
		Kpcb = 5.000000e-04;
		Kpc_max = 2.332400e-01;
		Kpc_half = 2.000000e+01;
		i_pCa_max = 1.000000e+00;
		Km_pCa = 5.000000e-01;
		k_NaCa = 2.928000e+02;
		K_mNa = 8.750000e+04;
		Nao = 1.400000e+05;
		K_mCa = 1.380000e+03;
		Cao = 1.800000e+03;
		k_sat = 1.000000e-01;
		eta = 3.500000e-01;
		R = 8.314000e+00;
		T = 2.980000e+02;
		g_Cab = 3.670000e-04;
//		g_Na = 1.300000e+01;
		Ko = 5.400000e+03;
		g_Nab = 2.600000e-03;
/*		g_Kto_f = 4.067000e-01;
		g_Kto_s = 0.000000e+00;
		g_Ks = 5.750000e-03;
		g_Kur = 1.600000e-01;
		g_Kss = 5.000000e-02;
		g_Kr = 7.800000e-02;*/
		kf = 2.376100e-02;
		kb = 3.677800e-02;
		i_NaK_max = 8.800000e-01;
		Km_Nai = 2.100000e+04;
		Km_Ko = 1.500000e+03;
		g_ClCa = 1.000000e+01;
		Km_Cl = 1.000000e+01;
		E_Cl = -4.000000e+01;
		IstimAmplitude = -8.000000e+01;
		IstimStart = 0.000000e+00;
		IstimEnd = 1.000000e+04;
		IstimPeriod = 1.000000e+03;
		IstimPulseDuration = 5.000000e-01;
                g_Kto_f = 0.0; // Inicializando com 0 para setar depois
                g_Kto_s = 0.0;
                g_Kur = 0.0;
                g_Kss = 0.0;
                g_Kr = 0.0;
                g_Ks = 0.0;
                g_Na = 0.0;
                g_CaL = 0.0;
		v_step = -8.26887e+01;
                protocolo = 0;
		dtime = 0.0;
                time_vec__ = NULL;

		dependent_variable__ = N_VNew_Serial(41);
		if(check_flag((void *)dependent_variable__, "N_VNew_Serial", 0))
			exit(1);
		abstol__ = N_VNew_Serial(41);
		if (check_flag((void *)abstol__, "N_VNew_Serial", 0))
			exit(1);

		depvar__ = (double*)malloc(sizeof(double)*41);
		if(depvar__ == NULL){
			fprintf(stderr, "ERROR Cannot allocate memory for depvar__\n");
			exit(0);
		}

                //condicoes iniciais da Dorothy
		V = NULL;
		//NV_Ith_S(dependent_variable__, 0) = V_ini_ = -8.242020e+01;	NV_Ith_S(abstol__, 0) = val_abstol__;
                NV_Ith_S(dependent_variable__, 0) = V_ini_ = -8.26887e+01;     NV_Ith_S(abstol__, 0) = val_abstol__;
		Cai = NULL;
		//NV_Ith_S(dependent_variable__, 1) = Cai_ini_ = 1.150010e-01;	NV_Ith_S(abstol__, 1) = val_abstol__;
                NV_Ith_S(dependent_variable__, 1) = Cai_ini_ = 1.06797e-01;    NV_Ith_S(abstol__, 1) = val_abstol__;
		Cass = NULL;
		//NV_Ith_S(dependent_variable__, 2) = Cass_ini_ = 1.150010e-01;	NV_Ith_S(abstol__, 2) = val_abstol__;
                NV_Ith_S(dependent_variable__, 2) = Cass_ini_ = 1.06797e-01;   NV_Ith_S(abstol__, 2) = val_abstol__;
		CaJSR = NULL;
		//NV_Ith_S(dependent_variable__, 3) = CaJSR_ini_ = 1.299500e+03;	NV_Ith_S(abstol__, 3) = val_abstol__;
                NV_Ith_S(dependent_variable__, 3) = CaJSR_ini_ = 1.03317e+03;  NV_Ith_S(abstol__, 3) = val_abstol__;
		CaNSR = NULL;
		//NV_Ith_S(dependent_variable__, 4) = CaNSR_ini_ = 1.299500e+03;	NV_Ith_S(abstol__, 4) = val_abstol__;
                NV_Ith_S(dependent_variable__, 4) = CaNSR_ini_ = 1.03469e+03;  NV_Ith_S(abstol__, 4) = val_abstol__;
		P_RyR = NULL;
		//NV_Ith_S(dependent_variable__, 5) = P_RyR_ini_ = 0.000000e+00;	NV_Ith_S(abstol__, 5) = val_abstol__;
                NV_Ith_S(dependent_variable__, 5) = P_RyR_ini_ = 4.18986e-22;  NV_Ith_S(abstol__, 5) = val_abstol__;
		LTRPN_Ca = NULL;
		//NV_Ith_S(dependent_variable__, 6) = LTRPN_Ca_ini_ = 1.126840e+01;	NV_Ith_S(abstol__, 6) = val_abstol__;
                NV_Ith_S(dependent_variable__, 6) = LTRPN_Ca_ini_ = 1.05919e+01;       NV_Ith_S(abstol__, 6) = val_abstol__;
		HTRPN_Ca = NULL;
		//NV_Ith_S(dependent_variable__, 7) = HTRPN_Ca_ini_ = 1.252900e+02;	NV_Ith_S(abstol__, 7) = val_abstol__;
                NV_Ith_S(dependent_variable__, 7) = HTRPN_Ca_ini_ = 1.2568e+02;       NV_Ith_S(abstol__, 7) = val_abstol__;
		P_O1 = NULL;
		//NV_Ith_S(dependent_variable__, 8) = P_O1_ini_ = 1.491020e-05;	NV_Ith_S(abstol__, 8) = val_abstol__;
                NV_Ith_S(dependent_variable__, 8) = P_O1_ini_ = 4.57256e-04; NV_Ith_S(abstol__, 8) = val_abstol__;
		P_O2 = NULL;
		//NV_Ith_S(dependent_variable__, 9) = P_O2_ini_ = 9.517260e-11;	NV_Ith_S(abstol__, 9) = val_abstol__;
                NV_Ith_S(dependent_variable__, 9) = P_O2_ini_ = 2.33937e-09;    NV_Ith_S(abstol__, 9) = val_abstol__;
		P_C2 = NULL;
		//NV_Ith_S(dependent_variable__, 10) = P_C2_ini_ = 1.677400e-04;	NV_Ith_S(abstol__, 10) = val_abstol__;
                NV_Ith_S(dependent_variable__, 10) = P_C2_ini_ = 4.45287e-02;   NV_Ith_S(abstol__, 10) = val_abstol__;
		O = NULL;
		//NV_Ith_S(dependent_variable__, 11) = O_ini_ = 9.303080e-19;	NV_Ith_S(abstol__, 11) = val_abstol__;
                NV_Ith_S(dependent_variable__, 11) = O_ini_ = 4.61989e-18;     NV_Ith_S(abstol__, 11) = val_abstol__;
		C2 = NULL;
		//NV_Ith_S(dependent_variable__, 12) = C2_ini_ = 1.242160e-04;	NV_Ith_S(abstol__, 12) = val_abstol__;
                NV_Ith_S(dependent_variable__, 12) = C2_ini_ = 1.18453e-04;    NV_Ith_S(abstol__, 12) = val_abstol__;
		C3 = NULL;
		//NV_Ith_S(dependent_variable__, 13) = C3_ini_ = 5.786790e-09;	NV_Ith_S(abstol__, 13) = val_abstol__;
                NV_Ith_S(dependent_variable__, 13) = C3_ini_ = 5.26227e-09;    NV_Ith_S(abstol__, 13) = val_abstol__;
		C4 = NULL;
		//NV_Ith_S(dependent_variable__, 14) = C4_ini_ = 1.198160e-13;	NV_Ith_S(abstol__, 14) = val_abstol__;
                NV_Ith_S(dependent_variable__, 14) = C4_ini_ = 1.04231e-13;    NV_Ith_S(abstol__, 14) = val_abstol__;
		I1 = NULL;
		//NV_Ith_S(dependent_variable__, 15) = I1_ini_ = 4.979230e-19;	NV_Ith_S(abstol__, 15) = val_abstol__;
                NV_Ith_S(dependent_variable__, 15) = I1_ini_ = 3.53862e-13;    NV_Ith_S(abstol__, 15) = val_abstol__;
		I2 = NULL;
		//NV_Ith_S(dependent_variable__, 16) = I2_ini_ = 3.458470e-14;	NV_Ith_S(abstol__, 16) = val_abstol__;
                NV_Ith_S(dependent_variable__, 16) = I2_ini_ = 3.24542e-14;    NV_Ith_S(abstol__, 16) = val_abstol__;
		I3 = NULL;
		//NV_Ith_S(dependent_variable__, 17) = I3_ini_ = 1.851060e-14;	NV_Ith_S(abstol__, 17) = val_abstol__;
                NV_Ith_S(dependent_variable__, 17) = I3_ini_ = 5.44951e-13;    NV_Ith_S(abstol__, 17) = val_abstol__;
		Nai = NULL;
		//NV_Ith_S(dependent_variable__, 18) = Nai_ini_ = 1.423710e+04;	NV_Ith_S(abstol__, 18) = val_abstol__;
                NV_Ith_S(dependent_variable__, 18) = Nai_ini_ = 1.59031e+04;   NV_Ith_S(abstol__, 18) = val_abstol__;
		C_Na2 = NULL;
		//NV_Ith_S(dependent_variable__, 19) = C_Na2_ini_ = 2.075200e-02;	NV_Ith_S(abstol__, 19) = val_abstol__;
                NV_Ith_S(dependent_variable__, 19) = C_Na2_ini_ = 2.04397e-02; NV_Ith_S(abstol__, 19) = val_abstol__;
		C_Na1 = NULL;
		//NV_Ith_S(dependent_variable__, 20) = C_Na1_ini_ = 2.791320e-04;	NV_Ith_S(abstol__, 20) = val_abstol__;
                NV_Ith_S(dependent_variable__, 20) = C_Na1_ini_ = 2.6658e-04; NV_Ith_S(abstol__, 20) = val_abstol__;
		O_Na = NULL;
		//NV_Ith_S(dependent_variable__, 21) = O_Na_ini_ = 7.134830e-07;	NV_Ith_S(abstol__, 21) = val_abstol__;
                NV_Ith_S(dependent_variable__, 21) = O_Na_ini_ = 6.57633e-07;  NV_Ith_S(abstol__, 21) = val_abstol__;
		IF_Na = NULL;
		//NV_Ith_S(dependent_variable__, 22) = IF_Na_ini_ = 1.531760e-04;	NV_Ith_S(abstol__, 22) = val_abstol__;
                NV_Ith_S(dependent_variable__, 22) = IF_Na_ini_ = 1.4115e-04; NV_Ith_S(abstol__, 22) = val_abstol__;
		I1_Na = NULL;
		//NV_Ith_S(dependent_variable__, 23) = I1_Na_ini_ = 6.733450e-07;	NV_Ith_S(abstol__, 23) = val_abstol__;
                NV_Ith_S(dependent_variable__, 23) = I1_Na_ini_ = 6.1422e-07; NV_Ith_S(abstol__, 23) = val_abstol__;
		I2_Na = NULL;
		//NV_Ith_S(dependent_variable__, 24) = I2_Na_ini_ = 1.557870e-09;	NV_Ith_S(abstol__, 24) = val_abstol__;
                NV_Ith_S(dependent_variable__, 24) = I2_Na_ini_ = 1.21421e-06; NV_Ith_S(abstol__, 24) = val_abstol__;
		IC_Na2 = NULL;
		//NV_Ith_S(dependent_variable__, 25) = IC_Na2_ini_ = 1.138790e-02;	NV_Ith_S(abstol__, 25) = val_abstol__;
                NV_Ith_S(dependent_variable__, 25) = IC_Na2_ini_ = 1.08225e-02;        NV_Ith_S(abstol__, 25) = val_abstol__;
		IC_Na3 = NULL;
		//NV_Ith_S(dependent_variable__, 26) = IC_Na3_ini_ = 3.427800e-01;	NV_Ith_S(abstol__, 26) = val_abstol__;
                NV_Ith_S(dependent_variable__, 26) = IC_Na3_ini_ = 3.35222e-01;        NV_Ith_S(abstol__, 26) = val_abstol__;
		Ki = NULL;
		//NV_Ith_S(dependent_variable__, 27) = Ki_ini_ = 1.437200e+05;	NV_Ith_S(abstol__, 27) = val_abstol__;
                NV_Ith_S(dependent_variable__, 27) = Ki_ini_ = 1.46344e+05;    NV_Ith_S(abstol__, 27) = val_abstol__;
		ato_f = NULL;
		NV_Ith_S(dependent_variable__, 28) = ato_f_ini_ = 2.655630e-03;	NV_Ith_S(abstol__, 28) = val_abstol__;
		ito_f = NULL;
		NV_Ith_S(dependent_variable__, 29) = ito_f_ini_ = 9.999770e-01;	NV_Ith_S(abstol__, 29) = val_abstol__;
		ato_s = NULL;
		NV_Ith_S(dependent_variable__, 30) = ato_s_ini_ = 4.170690e-04;	NV_Ith_S(abstol__, 30) = val_abstol__;
		ito_s = NULL;
		NV_Ith_S(dependent_variable__, 31) = ito_s_ini_ = 9.985430e-01;	NV_Ith_S(abstol__, 31) = val_abstol__;
		nKs = NULL;
		//NV_Ith_S(dependent_variable__, 32) = nKs_ini_ = 2.627530e-04;	NV_Ith_S(abstol__, 32) = val_abstol__;
                NV_Ith_S(dependent_variable__, 32) = nKs_ini_ = 4.20518e-04;   NV_Ith_S(abstol__, 32) = val_abstol__;
		aur = NULL;
		//NV_Ith_S(dependent_variable__, 33) = aur_ini_ = 4.170690e-04;	NV_Ith_S(abstol__, 33) = val_abstol__;
                NV_Ith_S(dependent_variable__, 33) = aur_ini_ = 4.02707e-04;   NV_Ith_S(abstol__, 33) = val_abstol__;
		iur = NULL;
		//NV_Ith_S(dependent_variable__, 34) = iur_ini_ = 9.985430e-01;	NV_Ith_S(abstol__, 34) = val_abstol__;
                NV_Ith_S(dependent_variable__, 34) = iur_ini_ = 9.97252e-01;   NV_Ith_S(abstol__, 34) = val_abstol__;
		aKss = NULL;
		//NV_Ith_S(dependent_variable__, 35) = aKss_ini_ = 4.170690e-04;	NV_Ith_S(abstol__, 35) = val_abstol__;
                NV_Ith_S(dependent_variable__, 35) = aKss_ini_ = 6.09841e-01;  NV_Ith_S(abstol__, 35) = val_abstol__;
		iKss = NULL;
		NV_Ith_S(dependent_variable__, 36) = iKss_ini_ = 1.000000e+00;	NV_Ith_S(abstol__, 36) = val_abstol__;
		C_K2 = NULL;
		//NV_Ith_S(dependent_variable__, 37) = C_K2_ini_ = 6.412290e-04;	NV_Ith_S(abstol__, 37) = val_abstol__;
                NV_Ith_S(dependent_variable__, 37) = C_K2_ini_ = 6.29406e-04;  NV_Ith_S(abstol__, 37) = val_abstol__;
		C_K1 = NULL;
		//NV_Ith_S(dependent_variable__, 38) = C_K1_ini_ = 9.925130e-04;	NV_Ith_S(abstol__, 38) = val_abstol__;
                NV_Ith_S(dependent_variable__, 38) = C_K1_ini_ = 9.72787e-04;  NV_Ith_S(abstol__, 38) = val_abstol__;
		O_K = NULL;
		//NV_Ith_S(dependent_variable__, 39) = O_K_ini_ = 1.752980e-04;	NV_Ith_S(abstol__, 39) = val_abstol__;
                NV_Ith_S(dependent_variable__, 39) = O_K_ini_ = 1.83473e-04;   NV_Ith_S(abstol__, 39) = val_abstol__;
		I_K = NULL;
		//NV_Ith_S(dependent_variable__, 40) = I_K_ini_ = 3.191290e-05;	NV_Ith_S(abstol__, 40) = val_abstol__;
                NV_Ith_S(dependent_variable__, 40) = I_K_ini_ = 3.29636e-05;   NV_Ith_S(abstol__, 40) = val_abstol__;
                istim = NULL;
                I_Ktof = NULL;
                I_Kur = NULL;
                I_Kss = NULL;
                I_Ks = NULL;
                I_Kr = NULL;
                I_K1 = NULL;
                I_CaL = NULL;
                I_Cab = NULL;
                I_Na = NULL;
                I_NaK = NULL;
                I_NaCa = NULL;
                I_Nab = NULL;
                I_pCa = NULL;
                I_ClCa = NULL;

		reltol__ = val_reltol__;
		it_countx = 0;
		int nonlineariteration__;
		if(method__ == CV_BDF) nonlineariteration__ = CV_NEWTON;
		else nonlineariteration__ = CV_FUNCTIONAL;
		cvode_mem_cvode__ = CVodeCreate(method__, nonlineariteration__);
		if (check_flag((void *)cvode_mem_cvode__, "CVodeCreate", 0))
			exit(1);

                flag__ = CVodeSetMaxStep(cvode_mem_cvode__, dtime);
                if (check_flag(&flag__, "MaxStep", 1))
                        exit(1);

		flag__ = CVodeMalloc(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SV, reltol__, abstol__);
		if (check_flag(&flag__, "CVodeMalloc", 1))
			exit(1);

/*		flag__ = CVDense(cvode_mem_cvode__, 41);
		if (check_flag(&flag__, "CVDense", 1))
			exit(1);*/

                flag__ = CVDiag(cvode_mem_cvode__);
                if (check_flag(&flag__, "CVDiag", 1))
                        exit(1);

		CVodeSetFdata(cvode_mem_cvode__, (void*)this);
	}
	Solveode::~Solveode()
	{
		if(depvar__ != NULL) free(depvar__);
		if(V != NULL) free(V);
		if(Cai != NULL) free(Cai);
		if(Cass != NULL) free(Cass);
		if(CaJSR != NULL) free(CaJSR);
		if(CaNSR != NULL) free(CaNSR);
		if(P_RyR != NULL) free(P_RyR);
		if(LTRPN_Ca != NULL) free(LTRPN_Ca);
		if(HTRPN_Ca != NULL) free(HTRPN_Ca);
		if(P_O1 != NULL) free(P_O1);
		if(P_O2 != NULL) free(P_O2);
		if(P_C2 != NULL) free(P_C2);
		if(O != NULL) free(O);
		if(C2 != NULL) free(C2);
		if(C3 != NULL) free(C3);
		if(C4 != NULL) free(C4);
		if(I1 != NULL) free(I1);
		if(I2 != NULL) free(I2);
		if(I3 != NULL) free(I3);
		if(Nai != NULL) free(Nai);
		if(C_Na2 != NULL) free(C_Na2);
		if(C_Na1 != NULL) free(C_Na1);
		if(O_Na != NULL) free(O_Na);
		if(IF_Na != NULL) free(IF_Na);
		if(I1_Na != NULL) free(I1_Na);
		if(I2_Na != NULL) free(I2_Na);
		if(IC_Na2 != NULL) free(IC_Na2);
		if(IC_Na3 != NULL) free(IC_Na3);
		if(Ki != NULL) free(Ki);
		if(ato_f != NULL) free(ato_f);
		if(ito_f != NULL) free(ito_f);
		if(ato_s != NULL) free(ato_s);
		if(ito_s != NULL) free(ito_s);
		if(nKs != NULL) free(nKs);
		if(aur != NULL) free(aur);
		if(iur != NULL) free(iur);
		if(aKss != NULL) free(aKss);
		if(iKss != NULL) free(iKss);
		if(C_K2 != NULL) free(C_K2);
		if(C_K1 != NULL) free(C_K1);
		if(O_K != NULL) free(O_K);
		if(I_K != NULL) free(I_K);

                if(istim != NULL) free(istim);
                if(I_Ktof != NULL) free(I_Ktof);
                if(I_Kur != NULL) free(I_Kur);
                if(I_Kss != NULL) free(I_Kss);
                if(I_Ks != NULL) free(I_Ks);
                if(I_Kr != NULL) free(I_Kr);
                if(I_K1 != NULL) free(I_K1);
                if(I_CaL != NULL) free(I_CaL);
                if(I_Cab != NULL) free(I_Cab);
                if(I_Na != NULL) free(I_Na);
                if(I_NaK != NULL) free(I_NaK);
                if(I_NaCa != NULL) free(I_NaCa);
                if(I_Nab != NULL) free(I_Nab);
                if(I_pCa != NULL) free(I_pCa);
                if(I_ClCa != NULL) free(I_ClCa);

		N_VDestroy_Serial(dependent_variable__);
		N_VDestroy_Serial(abstol__);
		CVodeFree(&cvode_mem_cvode__);
	}

	int Solveode::setVariables(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		NV_Ith_S(dependent_variable__, 0) = V_old_ = V_ini_ = value_new;    break;
		case 1:		NV_Ith_S(dependent_variable__, 1) = Cai_old_ = Cai_ini_ = value_new;    break;
		case 2:		NV_Ith_S(dependent_variable__, 2) = Cass_old_ = Cass_ini_ = value_new;    break;
		case 3:		NV_Ith_S(dependent_variable__, 3) = CaJSR_old_ = CaJSR_ini_ = value_new;    break;
		case 4:		NV_Ith_S(dependent_variable__, 4) = CaNSR_old_ = CaNSR_ini_ = value_new;    break;
		case 5:		NV_Ith_S(dependent_variable__, 5) = P_RyR_old_ = P_RyR_ini_ = value_new;    break;
		case 6:		NV_Ith_S(dependent_variable__, 6) = LTRPN_Ca_old_ = LTRPN_Ca_ini_ = value_new;    break;
		case 7:		NV_Ith_S(dependent_variable__, 7) = HTRPN_Ca_old_ = HTRPN_Ca_ini_ = value_new;    break;
		case 8:		NV_Ith_S(dependent_variable__, 8) = P_O1_old_ = P_O1_ini_ = value_new;    break;
		case 9:		NV_Ith_S(dependent_variable__, 9) = P_O2_old_ = P_O2_ini_ = value_new;    break;
		case 10:	NV_Ith_S(dependent_variable__, 10) = P_C2_old_ = P_C2_ini_ = value_new;    break;
		case 11:	NV_Ith_S(dependent_variable__, 11) = O_old_ = O_ini_ = value_new;    break;
		case 12:	NV_Ith_S(dependent_variable__, 12) = C2_old_ = C2_ini_ = value_new;    break;
		case 13:	NV_Ith_S(dependent_variable__, 13) = C3_old_ = C3_ini_ = value_new;    break;
		case 14:	NV_Ith_S(dependent_variable__, 14) = C4_old_ = C4_ini_ = value_new;    break;
		case 15:	NV_Ith_S(dependent_variable__, 15) = I1_old_ = I1_ini_ = value_new;    break;
		case 16:	NV_Ith_S(dependent_variable__, 16) = I2_old_ = I2_ini_ = value_new;    break;
		case 17:	NV_Ith_S(dependent_variable__, 17) = I3_old_ = I3_ini_ = value_new;    break;
		case 18:	NV_Ith_S(dependent_variable__, 18) = Nai_old_ = Nai_ini_ = value_new;    break;
		case 19:	NV_Ith_S(dependent_variable__, 19) = C_Na2_old_ = C_Na2_ini_ = value_new;    break;
		case 20:	NV_Ith_S(dependent_variable__, 20) = C_Na1_old_ = C_Na1_ini_ = value_new;    break;
		case 21:	NV_Ith_S(dependent_variable__, 21) = O_Na_old_ = O_Na_ini_ = value_new;    break;
		case 22:	NV_Ith_S(dependent_variable__, 22) = IF_Na_old_ = IF_Na_ini_ = value_new;    break;
		case 23:	NV_Ith_S(dependent_variable__, 23) = I1_Na_old_ = I1_Na_ini_ = value_new;    break;
		case 24:	NV_Ith_S(dependent_variable__, 24) = I2_Na_old_ = I2_Na_ini_ = value_new;    break;
		case 25:	NV_Ith_S(dependent_variable__, 25) = IC_Na2_old_ = IC_Na2_ini_ = value_new;    break;
		case 26:	NV_Ith_S(dependent_variable__, 26) = IC_Na3_old_ = IC_Na3_ini_ = value_new;    break;
		case 27:	NV_Ith_S(dependent_variable__, 27) = Ki_old_ = Ki_ini_ = value_new;    break;
		case 28:	NV_Ith_S(dependent_variable__, 28) = ato_f_old_ = ato_f_ini_ = value_new;    break;
		case 29:	NV_Ith_S(dependent_variable__, 29) = ito_f_old_ = ito_f_ini_ = value_new;    break;
		case 30:	NV_Ith_S(dependent_variable__, 30) = ato_s_old_ = ato_s_ini_ = value_new;    break;
		case 31:	NV_Ith_S(dependent_variable__, 31) = ito_s_old_ = ito_s_ini_ = value_new;    break;
		case 32:	NV_Ith_S(dependent_variable__, 32) = nKs_old_ = nKs_ini_ = value_new;    break;
		case 33:	NV_Ith_S(dependent_variable__, 33) = aur_old_ = aur_ini_ = value_new;    break;
		case 34:	NV_Ith_S(dependent_variable__, 34) = iur_old_ = iur_ini_ = value_new;    break;
		case 35:	NV_Ith_S(dependent_variable__, 35) = aKss_old_ = aKss_ini_ = value_new;    break;
		case 36:	NV_Ith_S(dependent_variable__, 36) = iKss_old_ = iKss_ini_ = value_new;    break;
		case 37:	NV_Ith_S(dependent_variable__, 37) = C_K2_old_ = C_K2_ini_ = value_new;    break;
		case 38:	NV_Ith_S(dependent_variable__, 38) = C_K1_old_ = C_K1_ini_ = value_new;    break;
		case 39:	NV_Ith_S(dependent_variable__, 39) = O_K_old_ = O_K_ini_ = value_new;    break;
		case 40:	NV_Ith_S(dependent_variable__, 40) = I_K_old_ = I_K_ini_ = value_new;    break;
		default:	return 1;    break;
		}
		return 0;
	}

	int Solveode::setParameters(int indVariable, double value_new)
	{
		switch(indVariable)
		{
		case 0:		time = value_new;   break;
		case 1:		i_stim = value_new;   break;
		case 2:		Cm = value_new;   break;
		case 3:		Acap = value_new;   break;
		case 4:		Vmyo = value_new;   break;
		case 5:		F = value_new;   break;
		case 6:		VJSR = value_new;   break;
		case 7:		Vss = value_new;   break;
		case 8:		VNSR = value_new;   break;
		case 9:		CMDN_tot = value_new;   break;
		case 10:		Km_CMDN = value_new;   break;
		case 11:		CSQN_tot = value_new;   break;
		case 12:		Km_CSQN = value_new;   break;
		case 13:		v1 = value_new;   break;
		case 14:		tau_tr = value_new;   break;
		case 15:		tau_xfer = value_new;   break;
		case 16:		v2 = value_new;   break;
		case 17:		v3 = value_new;   break;
		case 18:		Km_up = value_new;   break;
		case 19:		k_plus_htrpn = value_new;   break;
		case 20:		HTRPN_tot = value_new;   break;
		case 21:		k_plus_ltrpn = value_new;   break;
		case 22:		LTRPN_tot = value_new;   break;
		case 23:		k_minus_htrpn = value_new;   break;
		case 24:		k_minus_ltrpn = value_new;   break;
		case 25:		i_CaL_max = value_new;   break;
		case 26:		k_plus_a = value_new;   break;
		case 27:		n = value_new;   break;
		case 28:		k_minus_b = value_new;   break;
		case 29:		k_minus_c = value_new;   break;
		case 30:		k_minus_a = value_new;   break;
		case 31:		k_plus_b = value_new;   break;
		case 32:		m = value_new;   break;
		case 33:		k_plus_c = value_new;   break;
		case 34:		g_CaL = value_new;   break;
		case 35:		E_CaL = value_new;   break;
		case 36:		Kpcb = value_new;   break;
		case 37:		Kpc_max = value_new;   break;
		case 38:		Kpc_half = value_new;   break;
		case 39:		i_pCa_max = value_new;   break;
		case 40:		Km_pCa = value_new;   break;
		case 41:		k_NaCa = value_new;   break;
		case 42:		K_mNa = value_new;   break;
		case 43:		Nao = value_new;   break;
		case 44:		K_mCa = value_new;   break;
		case 45:		Cao = value_new;   break;
		case 46:		k_sat = value_new;   break;
		case 47:		eta = value_new;   break;
		case 48:		R = value_new;   break;
		case 49:		T = value_new;   break;
		case 50:		g_Cab = value_new;   break;
		case 51:		g_Na = value_new;   break;
		case 52:		Ko = value_new;   break;
		case 53:		g_Nab = value_new;   break;
		case 54:		g_Kto_f = value_new;   break;
		case 55:		g_Kto_s = value_new;   break;
		case 56:		g_Ks = value_new;   break;
		case 57:		g_Kur = value_new;   break;
		case 58:		g_Kss = value_new;   break;
		case 59:		g_Kr = value_new;   break;
		case 60:		kf = value_new;   break;
		case 61:		kb = value_new;   break;
		case 62:		i_NaK_max = value_new;   break;
		case 63:		Km_Nai = value_new;   break;
		case 64:		Km_Ko = value_new;   break;
		case 65:		g_ClCa = value_new;   break;
		case 66:		Km_Cl = value_new;   break;
		case 67:		E_Cl = value_new;   break;
		case 68:		IstimAmplitude = value_new;   break;
		case 69:		IstimStart = value_new;   break;
		case 70:		IstimEnd = value_new;   break;
		case 71:		IstimPeriod = value_new;   break;
		case 72:		IstimPulseDuration = value_new;   break;
                case 73:                v_step = value_new;   break;
                case 74:                protocolo = value_new;   break;
		default:	return 1;    break;
		}
		return 0;
	}

	int Solveode::setFreeVariable(double value_new)
	{
		dtime = value_new;
	}

	//Get Methods

	double Solveode::getVariables(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V_ini_;    break;
		case 1:		return Cai_ini_;    break;
		case 2:		return Cass_ini_;    break;
		case 3:		return CaJSR_ini_;    break;
		case 4:		return CaNSR_ini_;    break;
		case 5:		return P_RyR_ini_;    break;
		case 6:		return LTRPN_Ca_ini_;    break;
		case 7:		return HTRPN_Ca_ini_;    break;
		case 8:		return P_O1_ini_;    break;
		case 9:		return P_O2_ini_;    break;
		case 10:		return P_C2_ini_;    break;
		case 11:		return O_ini_;    break;
		case 12:		return C2_ini_;    break;
		case 13:		return C3_ini_;    break;
		case 14:		return C4_ini_;    break;
		case 15:		return I1_ini_;    break;
		case 16:		return I2_ini_;    break;
		case 17:		return I3_ini_;    break;
		case 18:		return Nai_ini_;    break;
		case 19:		return C_Na2_ini_;    break;
		case 20:		return C_Na1_ini_;    break;
		case 21:		return O_Na_ini_;    break;
		case 22:		return IF_Na_ini_;    break;
		case 23:		return I1_Na_ini_;    break;
		case 24:		return I2_Na_ini_;    break;
		case 25:		return IC_Na2_ini_;    break;
		case 26:		return IC_Na3_ini_;    break;
		case 27:		return Ki_ini_;    break;
		case 28:		return ato_f_ini_;    break;
		case 29:		return ito_f_ini_;    break;
		case 30:		return ato_s_ini_;    break;
		case 31:		return ito_s_ini_;    break;
		case 32:		return nKs_ini_;    break;
		case 33:		return aur_ini_;    break;
		case 34:		return iur_ini_;    break;
		case 35:		return aKss_ini_;    break;
		case 36:		return iKss_ini_;    break;
		case 37:		return C_K2_ini_;    break;
		case 38:		return C_K1_ini_;    break;
		case 39:		return O_K_ini_;    break;
		case 40:		return I_K_ini_;    break;
		default:	return 1;    break;
		}
	}

	double Solveode::getParameters(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return time;    break;
		case 1:		return i_stim;    break;
		case 2:		return Cm;    break;
		case 3:		return Acap;    break;
		case 4:		return Vmyo;    break;
		case 5:		return F;    break;
		case 6:		return VJSR;    break;
		case 7:		return Vss;    break;
		case 8:		return VNSR;    break;
		case 9:		return CMDN_tot;    break;
		case 10:		return Km_CMDN;    break;
		case 11:		return CSQN_tot;    break;
		case 12:		return Km_CSQN;    break;
		case 13:		return v1;    break;
		case 14:		return tau_tr;    break;
		case 15:		return tau_xfer;    break;
		case 16:		return v2;    break;
		case 17:		return v3;    break;
		case 18:		return Km_up;    break;
		case 19:		return k_plus_htrpn;    break;
		case 20:		return HTRPN_tot;    break;
		case 21:		return k_plus_ltrpn;    break;
		case 22:		return LTRPN_tot;    break;
		case 23:		return k_minus_htrpn;    break;
		case 24:		return k_minus_ltrpn;    break;
		case 25:		return i_CaL_max;    break;
		case 26:		return k_plus_a;    break;
		case 27:		return n;    break;
		case 28:		return k_minus_b;    break;
		case 29:		return k_minus_c;    break;
		case 30:		return k_minus_a;    break;
		case 31:		return k_plus_b;    break;
		case 32:		return m;    break;
		case 33:		return k_plus_c;    break;
		case 34:		return g_CaL;    break;
		case 35:		return E_CaL;    break;
		case 36:		return Kpcb;    break;
		case 37:		return Kpc_max;    break;
		case 38:		return Kpc_half;    break;
		case 39:		return i_pCa_max;    break;
		case 40:		return Km_pCa;    break;
		case 41:		return k_NaCa;    break;
		case 42:		return K_mNa;    break;
		case 43:		return Nao;    break;
		case 44:		return K_mCa;    break;
		case 45:		return Cao;    break;
		case 46:		return k_sat;    break;
		case 47:		return eta;    break;
		case 48:		return R;    break;
		case 49:		return T;    break;
		case 50:		return g_Cab;    break;
		case 51:		return g_Na;    break;
		case 52:		return Ko;    break;
		case 53:		return g_Nab;    break;
		case 54:		return g_Kto_f;    break;
		case 55:		return g_Kto_s;    break;
		case 56:		return g_Ks;    break;
		case 57:		return g_Kur;    break;
		case 58:		return g_Kss;    break;
		case 59:		return g_Kr;    break;
		case 60:		return kf;    break;
		case 61:		return kb;    break;
		case 62:		return i_NaK_max;    break;
		case 63:		return Km_Nai;    break;
		case 64:		return Km_Ko;    break;
		case 65:		return g_ClCa;    break;
		case 66:		return Km_Cl;    break;
		case 67:		return E_Cl;    break;
		case 68:		return IstimAmplitude;    break;
		case 69:		return IstimStart;    break;
		case 70:		return IstimEnd;    break;
		case 71:		return IstimPeriod;    break;
		case 72:		return IstimPulseDuration;    break;
                case 73:                return v_step;    break;
                case 74:                return protocolo;    break;
		default:	break;
		}
	}

	double Solveode::getFreeVariable()
	{
		return dtime;
	}

	//Get Methods - Variables

	Variables Solveode::get_Variables()
	{
		Variables v("|V#|Cai#|Cass#|CaJSR#|CaNSR#|P_RyR#|LTRPN_Ca#|HTRPN_Ca#|P_O1#|P_O2#|P_C2#|O#|C2#|C3#|C4#|I1#|I2#|I3#|Nai#|C_Na2#|C_Na1#|O_Na#|IF_Na#|I1_Na#|I2_Na#|IC_Na2#|IC_Na3#|Ki#|ato_f#|ito_f#|ato_s#|ito_s#|nKs#|aur#|iur#|aKss#|iKss#|C_K2#|C_K1#|O_K#|I_K#");
		return v;
	}
	Variables Solveode::get_Parameters()
	{
		Variables v("|time#|i_stim#|Cm#|Acap#|Vmyo#|F#|VJSR#|Vss#|VNSR#|CMDN_tot#|Km_CMDN#|CSQN_tot#|Km_CSQN#|v1#|tau_tr#|tau_xfer#|v2#|v3#|Km_up#|k_plus_htrpn#|HTRPN_tot#|k_plus_ltrpn#|LTRPN_tot#|k_minus_htrpn#|k_minus_ltrpn#|i_CaL_max#|k_plus_a#|n#|k_minus_b#|k_minus_c#|k_minus_a#|k_plus_b#|m#|k_plus_c#|g_CaL#|E_CaL#|Kpcb#|Kpc_max#|Kpc_half#|i_pCa_max#|Km_pCa#|k_NaCa#|K_mNa#|Nao#|K_mCa#|Cao#|k_sat#|eta#|R#|T#|g_Cab#|g_Na#|Ko#|g_Nab#|g_Kto_f#|g_Kto_s#|g_Ks#|g_Kur#|g_Kss#|g_Kr#|kf#|kb#|i_NaK_max#|Km_Nai#|Km_Ko#|g_ClCa#|Km_Cl#|E_Cl#|IstimAmplitude#|IstimStart#|IstimEnd#|IstimPeriod#|IstimPulseDuration#");
		return v;
	}
	Variables Solveode::get_FreeVariable()
	{
		Variables v("|time#");
		return v;
	}

	void Solveode::setParametersFromFile(char *filename)
	{
		FILE *file;
		if((file = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr,"ERROR - setParametersFromFile - Unable to open file %s\n", filename);
			exit(1);
		}
		double value;
		int k = 0;
		Variables v = get_Parameters();
		int s = v.getQuantity();
		for(;k<s;k++)
		{
			fscanf(file,"%lf", &value);
			setParameters(k, value);
		}
		fclose(file);
	}

	void Solveode::setVariablesFromFile(char *filename)
	{
		FILE *file;
		if((file = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr,"ERROR - setVariablesFromFile - Unable to open file %s\n", filename);
			exit(1);
		}
		double value;
		int k = 0;
		Variables v = get_Variables();
		int s = v.getQuantity();
		for(;k<s;k++)
		{
			fscanf(file,"%lf", &value);
			setVariables(k, value);
		}
		fclose(file);
	}

	void Solveode::setFreeVariableFromFile(char *filename)
	{
		FILE *file;
		if((file = fopen(filename, "r")) == NULL)
		{
			fprintf(stderr,"ERROR - setFreeVariableFromFile - Unable to open file %s\n", filename);
			exit(1);
		}
		double value;
		fscanf(file,"%lf", &value);
			setFreeVariable(value);
		fclose(file);
	}

	double* Solveode::solveDiff()
	{

		depvar__[0]= ((-(calc_i_CaL()+calc_i_pCa()+calc_i_NaCa()+calc_i_Cab()+calc_i_Na()+calc_i_Nab()+calc_i_NaK()+calc_i_Kto_f()+calc_i_Kto_s()+calc_i_K1()+calc_i_Ks()+calc_i_Kur()+calc_i_Kss()+calc_i_Kr()+calc_i_ClCa()+calc_Istim()))/Cm);	// 0
		depvar__[1]= (calc_Bi()*((calc_J_leak()+calc_J_xfer())-(calc_J_up()+calc_J_trpn()+((((calc_i_Cab()+calc_i_pCa())-(2.*calc_i_NaCa()))*Acap*Cm)/(2.*Vmyo*F)))));	// 1
		depvar__[2]= (calc_Bss()*(((calc_J_rel()*VJSR)/Vss)-(((calc_J_xfer()*Vmyo)/Vss)+((calc_i_CaL()*Acap*Cm)/(2.*Vss*F)))));	// 2
		depvar__[3]= (calc_BJSR()*(calc_J_tr()-calc_J_rel()));	// 3
		depvar__[4]= ((((calc_J_up()-calc_J_leak())*Vmyo)/VNSR)-((calc_J_tr()*VJSR)/VNSR));	// 4
		depvar__[5]= (((-0.04)*P_RyR_old_)-(((0.1*calc_i_CaL())/i_CaL_max)*exp(((-pow((V_old_-5.),2.))/648.))));	// 14
		depvar__[6]= ((k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_))-(k_minus_ltrpn*LTRPN_Ca_old_));	// 15
		depvar__[7]= ((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))-(k_minus_htrpn*HTRPN_Ca_old_));	// 16
		depvar__[8]= (((k_plus_a*pow(Cass_old_,n)*calc_P_C1())+(k_minus_b*P_O2_old_)+(k_minus_c*P_C2_old_))-((k_minus_a*P_O1_old_)+(k_plus_b*pow(Cass_old_,m)*P_O1_old_)+(k_plus_c*P_O1_old_)));	// 17
		depvar__[9]= ((k_plus_b*pow(Cass_old_,m)*P_O1_old_)-(k_minus_b*P_O2_old_));	// 19
		depvar__[10]= ((k_plus_c*P_O1_old_)-(k_minus_c*P_C2_old_));	// 20
		depvar__[11]= (((calc_alpha()*C4_old_)+(Kpcb*I1_old_)+(0.001*((calc_alpha()*I2_old_)-(calc_Kpcf()*O_old_))))-((4.*calc_beta()*O_old_)+(calc_gamma()*O_old_)));	// 22
		depvar__[12]= (((4.*calc_alpha()*calc_C1())+(2.*calc_beta()*C3_old_))-((calc_beta()*C2_old_)+(3.*calc_alpha()*C2_old_)));	// 24
		depvar__[13]= (((3.*calc_alpha()*C2_old_)+(3.*calc_beta()*C4_old_))-((2.*calc_beta()*C3_old_)+(2.*calc_alpha()*C3_old_)));	// 25
		depvar__[14]= (((2.*calc_alpha()*C3_old_)+(4.*calc_beta()*O_old_)+(0.01*((4.*Kpcb*calc_beta()*I1_old_)-(calc_alpha()*calc_gamma()*C4_old_)))+(0.002*((4.*calc_beta()*I2_old_)-(calc_Kpcf()*C4_old_)))+(4.*calc_beta()*Kpcb*I3_old_))-((3.*calc_beta()*C4_old_)+(calc_alpha()*C4_old_)+(calc_gamma()*calc_Kpcf()*C4_old_)));	// 26
		depvar__[15]= (((calc_gamma()*O_old_)+(0.001*((calc_alpha()*I3_old_)-(calc_Kpcf()*I1_old_)))+(0.01*((calc_alpha()*calc_gamma()*C4_old_)-(4.*calc_beta()*calc_Kpcf()*I1_old_))))-(Kpcb*I1_old_));	// 27
		depvar__[16]= (((0.001*((calc_Kpcf()*O_old_)-(calc_alpha()*I2_old_)))+(Kpcb*I3_old_)+(0.002*((calc_Kpcf()*C4_old_)-(4.*calc_beta()*I2_old_))))-(calc_gamma()*I2_old_));	// 28
		depvar__[17]= (((0.001*((calc_Kpcf()*I1_old_)-(calc_alpha()*I3_old_)))+(calc_gamma()*I2_old_)+(calc_gamma()*calc_Kpcf()*C4_old_))-((4.*calc_beta()*Kpcb*I3_old_)+(Kpcb*I3_old_)));	// 29
		depvar__[18]= (((-(calc_i_Na()+calc_i_Nab()+(3.*calc_i_NaK())+(3.*calc_i_NaCa())))*Acap*Cm)/(Vmyo*F));	// 38
		depvar__[19]= (((calc_alpha_Na11()*calc_C_Na3())+(calc_beta_Na12()*C_Na1_old_)+(calc_alpha_Na3()*IC_Na2_old_))-((calc_beta_Na11()*C_Na2_old_)+(calc_alpha_Na12()*C_Na2_old_)+(calc_beta_Na3()*C_Na2_old_)));	// 42
		depvar__[20]= (((calc_alpha_Na12()*C_Na2_old_)+(calc_beta_Na13()*O_Na_old_)+(calc_alpha_Na3()*IF_Na_old_))-((calc_beta_Na12()*C_Na1_old_)+(calc_alpha_Na13()*C_Na1_old_)+(calc_beta_Na3()*C_Na1_old_)));	// 43
		depvar__[21]= (((calc_alpha_Na13()*C_Na1_old_)+(calc_beta_Na2()*IF_Na_old_))-((calc_beta_Na13()*O_Na_old_)+(calc_alpha_Na2()*O_Na_old_)));	// 44
		depvar__[22]= (((calc_alpha_Na2()*O_Na_old_)+(calc_beta_Na3()*C_Na1_old_)+(calc_beta_Na4()*I1_Na_old_)+(calc_alpha_Na12()*IC_Na2_old_))-((calc_beta_Na2()*IF_Na_old_)+(calc_alpha_Na3()*IF_Na_old_)+(calc_alpha_Na4()*IF_Na_old_)+(calc_beta_Na12()*IF_Na_old_)));	// 45
		depvar__[23]= (((calc_alpha_Na4()*IF_Na_old_)+(calc_beta_Na5()*I2_Na_old_))-((calc_beta_Na4()*I1_Na_old_)+(calc_alpha_Na5()*I1_Na_old_)));	// 46
		depvar__[24]= ((calc_alpha_Na5()*I1_Na_old_)-(calc_beta_Na5()*I2_Na_old_));	// 47
		depvar__[25]= (((calc_alpha_Na11()*IC_Na3_old_)+(calc_beta_Na12()*IF_Na_old_)+(calc_beta_Na3()*IC_Na2_old_))-((calc_beta_Na11()*IC_Na2_old_)+(calc_alpha_Na12()*IC_Na2_old_)+(calc_alpha_Na3()*IC_Na2_old_)));	// 48
		depvar__[26]= (((calc_beta_Na11()*IC_Na2_old_)+(calc_beta_Na3()*calc_C_Na3()))-((calc_alpha_Na11()*IC_Na3_old_)+(calc_alpha_Na3()*IC_Na3_old_)));	// 49
		depvar__[27]= (((-((calc_Istim()+calc_i_Kto_f()+calc_i_Kto_s()+calc_i_K1()+calc_i_Ks()+calc_i_Kss()+calc_i_Kur()+calc_i_Kr())-(2.*calc_i_NaK())))*Acap*Cm)/(Vmyo*F));	// 65
		depvar__[28]= ((calc_alpha_a()*(1.-ato_f_old_))-(calc_beta_a()*ato_f_old_));	// 68
		depvar__[29]= ((calc_alpha_i()*(1.-ito_f_old_))-(calc_beta_i()*ito_f_old_));	// 69
		depvar__[30]= ((calc_ass()-ato_s_old_)/calc_tau_ta_s());	// 75
		depvar__[31]= ((calc_iss()-ito_s_old_)/calc_tau_ti_s());	// 76
		depvar__[32]= ((calc_alpha_n()*(1.-nKs_old_))-(calc_beta_n()*nKs_old_));	// 83
		depvar__[33]= ((calc_ass()-aur_old_)/calc_tau_aur());	// 87
		depvar__[34]= ((calc_iss()-iur_old_)/calc_tau_iur());	// 88
		depvar__[35]= ((calc_ass()-aKss_old_)/calc_tau_Kss());	// 92
		depvar__[36]= 0.;	// 93
		depvar__[37]= (((kf*C_K1_old_)+(calc_beta_a1()*O_K_old_))-((kb*C_K2_old_)+(calc_alpha_a1()*C_K2_old_)));	// 97
		depvar__[38]= (((calc_alpha_a0()*calc_C_K0())+(kb*C_K2_old_))-((calc_beta_a0()*C_K1_old_)+(kf*C_K1_old_)));	// 98
		depvar__[39]= (((calc_alpha_a1()*C_K2_old_)+(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current()*I_K_old_))-((calc_beta_a1()*O_K_old_)+(calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current()*O_K_old_)));	// 99
		depvar__[40]= ((calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current()*O_K_old_)-(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current()*I_K_old_));	// 100
		return depvar__;
	}


	void Solveode::solveCVODE(int firstcall__, int steps__)
	{
		int iout = 0;
		realtype tout = time_new+dtime;

		static int num_iterations_bak = 0;
		if(firstcall__){
			time_new = time;

			tout = time + dtime;
			if(steps__ <= 0)
				steps__ = 1;
			if(time_vec__ != NULL)free( time_vec__);
			time_vec__ = (double *)malloc(sizeof(double)*steps__);
			V_old_ = V_ini_;
			if(V != NULL)free( V);
			V = (double *)malloc(sizeof(double)*steps__);
			Cai_old_ = Cai_ini_;
			if(Cai != NULL)free( Cai);
			Cai = (double *)malloc(sizeof(double)*steps__);
			Cass_old_ = Cass_ini_;
			if(Cass != NULL)free( Cass);
			Cass = (double *)malloc(sizeof(double)*steps__);
			CaJSR_old_ = CaJSR_ini_;
			if(CaJSR != NULL)free( CaJSR);
			CaJSR = (double *)malloc(sizeof(double)*steps__);
			CaNSR_old_ = CaNSR_ini_;
			if(CaNSR != NULL)free( CaNSR);
			CaNSR = (double *)malloc(sizeof(double)*steps__);
			P_RyR_old_ = P_RyR_ini_;
			if(P_RyR != NULL)free( P_RyR);
			P_RyR = (double *)malloc(sizeof(double)*steps__);
			LTRPN_Ca_old_ = LTRPN_Ca_ini_;
			if(LTRPN_Ca != NULL)free( LTRPN_Ca);
			LTRPN_Ca = (double *)malloc(sizeof(double)*steps__);
			HTRPN_Ca_old_ = HTRPN_Ca_ini_;
			if(HTRPN_Ca != NULL)free( HTRPN_Ca);
			HTRPN_Ca = (double *)malloc(sizeof(double)*steps__);
			P_O1_old_ = P_O1_ini_;
			if(P_O1 != NULL)free( P_O1);
			P_O1 = (double *)malloc(sizeof(double)*steps__);
			P_O2_old_ = P_O2_ini_;
			if(P_O2 != NULL)free( P_O2);
			P_O2 = (double *)malloc(sizeof(double)*steps__);
			P_C2_old_ = P_C2_ini_;
			if(P_C2 != NULL)free( P_C2);
			P_C2 = (double *)malloc(sizeof(double)*steps__);
			O_old_ = O_ini_;
			if(O != NULL)free( O);
			O = (double *)malloc(sizeof(double)*steps__);
			C2_old_ = C2_ini_;
			if(C2 != NULL)free( C2);
			C2 = (double *)malloc(sizeof(double)*steps__);
			C3_old_ = C3_ini_;
			if(C3 != NULL)free( C3);
			C3 = (double *)malloc(sizeof(double)*steps__);
			C4_old_ = C4_ini_;
			if(C4 != NULL)free( C4);
			C4 = (double *)malloc(sizeof(double)*steps__);
			I1_old_ = I1_ini_;
			if(I1 != NULL)free( I1);
			I1 = (double *)malloc(sizeof(double)*steps__);
			I2_old_ = I2_ini_;
			if(I2 != NULL)free( I2);
			I2 = (double *)malloc(sizeof(double)*steps__);
			I3_old_ = I3_ini_;
			if(I3 != NULL)free( I3);
			I3 = (double *)malloc(sizeof(double)*steps__);
			Nai_old_ = Nai_ini_;
			if(Nai != NULL)free( Nai);
			Nai = (double *)malloc(sizeof(double)*steps__);
			C_Na2_old_ = C_Na2_ini_;
			if(C_Na2 != NULL)free( C_Na2);
			C_Na2 = (double *)malloc(sizeof(double)*steps__);
			C_Na1_old_ = C_Na1_ini_;
			if(C_Na1 != NULL)free( C_Na1);
			C_Na1 = (double *)malloc(sizeof(double)*steps__);
			O_Na_old_ = O_Na_ini_;
			if(O_Na != NULL)free( O_Na);
			O_Na = (double *)malloc(sizeof(double)*steps__);
			IF_Na_old_ = IF_Na_ini_;
			if(IF_Na != NULL)free( IF_Na);
			IF_Na = (double *)malloc(sizeof(double)*steps__);
			I1_Na_old_ = I1_Na_ini_;
			if(I1_Na != NULL)free( I1_Na);
			I1_Na = (double *)malloc(sizeof(double)*steps__);
			I2_Na_old_ = I2_Na_ini_;
			if(I2_Na != NULL)free( I2_Na);
			I2_Na = (double *)malloc(sizeof(double)*steps__);
			IC_Na2_old_ = IC_Na2_ini_;
			if(IC_Na2 != NULL)free( IC_Na2);
			IC_Na2 = (double *)malloc(sizeof(double)*steps__);
			IC_Na3_old_ = IC_Na3_ini_;
			if(IC_Na3 != NULL)free( IC_Na3);
			IC_Na3 = (double *)malloc(sizeof(double)*steps__);
			Ki_old_ = Ki_ini_;
			if(Ki != NULL)free( Ki);
			Ki = (double *)malloc(sizeof(double)*steps__);
			ato_f_old_ = ato_f_ini_;
			if(ato_f != NULL)free( ato_f);
			ato_f = (double *)malloc(sizeof(double)*steps__);
			ito_f_old_ = ito_f_ini_;
			if(ito_f != NULL)free( ito_f);
			ito_f = (double *)malloc(sizeof(double)*steps__);
			ato_s_old_ = ato_s_ini_;
			if(ato_s != NULL)free( ato_s);
			ato_s = (double *)malloc(sizeof(double)*steps__);
			ito_s_old_ = ito_s_ini_;
			if(ito_s != NULL)free( ito_s);
			ito_s = (double *)malloc(sizeof(double)*steps__);
			nKs_old_ = nKs_ini_;
			if(nKs != NULL)free( nKs);
			nKs = (double *)malloc(sizeof(double)*steps__);
			aur_old_ = aur_ini_;
			if(aur != NULL)free( aur);
			aur = (double *)malloc(sizeof(double)*steps__);
			iur_old_ = iur_ini_;
			if(iur != NULL)free( iur);
			iur = (double *)malloc(sizeof(double)*steps__);
			aKss_old_ = aKss_ini_;
			if(aKss != NULL)free( aKss);
			aKss = (double *)malloc(sizeof(double)*steps__);
			iKss_old_ = iKss_ini_;
			if(iKss != NULL)free( iKss);
			iKss = (double *)malloc(sizeof(double)*steps__);
			C_K2_old_ = C_K2_ini_;
			if(C_K2 != NULL)free( C_K2);
			C_K2 = (double *)malloc(sizeof(double)*steps__);
			C_K1_old_ = C_K1_ini_;
			if(C_K1 != NULL)free( C_K1);
			C_K1 = (double *)malloc(sizeof(double)*steps__);
			O_K_old_ = O_K_ini_;
			if(O_K != NULL)free( O_K);
			O_K = (double *)malloc(sizeof(double)*steps__);
			I_K_old_ = I_K_ini_;
			if(I_K != NULL)free( I_K);
			I_K = (double *)malloc(sizeof(double)*steps__);
                        if(istim != NULL)free(istim);
                        istim = (double *)malloc(sizeof(double)*steps__);

                        if(I_Ktof != NULL)free(I_Ktof);
                        I_Ktof = (double *)malloc(sizeof(double)*steps__);
                        if(I_Kur != NULL)free(I_Kur);
                        I_Kur = (double *)malloc(sizeof(double)*steps__);
                        if(I_Kss != NULL)free(I_Kss);
                        I_Kss = (double *)malloc(sizeof(double)*steps__);
                        if(I_Ks != NULL)free(I_Ks);
                        I_Ks = (double *)malloc(sizeof(double)*steps__);
                        if(I_Kr != NULL)free(I_Kr);
                        I_Kr = (double *)malloc(sizeof(double)*steps__);
                        if(I_K1 != NULL)free(I_K1);
                        I_K1 = (double *)malloc(sizeof(double)*steps__);
                        if(I_CaL != NULL)free(I_CaL);
                        I_CaL = (double *)malloc(sizeof(double)*steps__);
                        if(I_Cab != NULL)free(I_Cab);
                        I_Cab = (double *)malloc(sizeof(double)*steps__);
                        if(I_Na != NULL)free(I_Na);
                        I_Na = (double *)malloc(sizeof(double)*steps__);
                        if(I_NaK != NULL)free(I_NaK);
                        I_NaK = (double *)malloc(sizeof(double)*steps__);
                        if(I_NaCa != NULL)free(I_NaCa);
                        I_NaCa = (double *)malloc(sizeof(double)*steps__);
                        if(I_Nab != NULL)free(I_Nab);
                        I_Nab = (double *)malloc(sizeof(double)*steps__);
                        if(I_pCa != NULL)free(I_pCa);
                        I_pCa = (double *)malloc(sizeof(double)*steps__);
                        if(I_ClCa != NULL)free(I_ClCa);
                        I_ClCa = (double *)malloc(sizeof(double)*steps__);

			num_iterations_bak = steps__;
		}
		while(1){
                        flag__ = CVodeSetMaxStep(cvode_mem_cvode__, dtime);
                        flag__ = CVode(cvode_mem_cvode__, tout ,dependent_variable__, &time_new, CV_NORMAL);
                        flag__ = CVodeSetMaxStep(cvode_mem_cvode__, dtime);
			time_vec__[iout] = time_new;
			V_old_ = V[iout] = NV_Ith_S(dependent_variable__, 0);
			Cai_old_ = Cai[iout] = NV_Ith_S(dependent_variable__, 1);
			Cass_old_ = Cass[iout] = NV_Ith_S(dependent_variable__, 2);
			CaJSR_old_ = CaJSR[iout] = NV_Ith_S(dependent_variable__, 3);
			CaNSR_old_ = CaNSR[iout] = NV_Ith_S(dependent_variable__, 4);
			P_RyR_old_ = P_RyR[iout] = NV_Ith_S(dependent_variable__, 5);
			LTRPN_Ca_old_ = LTRPN_Ca[iout] = NV_Ith_S(dependent_variable__, 6);
			HTRPN_Ca_old_ = HTRPN_Ca[iout] = NV_Ith_S(dependent_variable__, 7);
			P_O1_old_ = P_O1[iout] = NV_Ith_S(dependent_variable__, 8);
			P_O2_old_ = P_O2[iout] = NV_Ith_S(dependent_variable__, 9);
			P_C2_old_ = P_C2[iout] = NV_Ith_S(dependent_variable__, 10);
			O_old_ = O[iout] = NV_Ith_S(dependent_variable__, 11);
			C2_old_ = C2[iout] = NV_Ith_S(dependent_variable__, 12);
			C3_old_ = C3[iout] = NV_Ith_S(dependent_variable__, 13);
			C4_old_ = C4[iout] = NV_Ith_S(dependent_variable__, 14);
			I1_old_ = I1[iout] = NV_Ith_S(dependent_variable__, 15);
			I2_old_ = I2[iout] = NV_Ith_S(dependent_variable__, 16);
			I3_old_ = I3[iout] = NV_Ith_S(dependent_variable__, 17);
			Nai_old_ = Nai[iout] = NV_Ith_S(dependent_variable__, 18);
			C_Na2_old_ = C_Na2[iout] = NV_Ith_S(dependent_variable__, 19);
			C_Na1_old_ = C_Na1[iout] = NV_Ith_S(dependent_variable__, 20);
			O_Na_old_ = O_Na[iout] = NV_Ith_S(dependent_variable__, 21);
			IF_Na_old_ = IF_Na[iout] = NV_Ith_S(dependent_variable__, 22);
			I1_Na_old_ = I1_Na[iout] = NV_Ith_S(dependent_variable__, 23);
			I2_Na_old_ = I2_Na[iout] = NV_Ith_S(dependent_variable__, 24);
			IC_Na2_old_ = IC_Na2[iout] = NV_Ith_S(dependent_variable__, 25);
			IC_Na3_old_ = IC_Na3[iout] = NV_Ith_S(dependent_variable__, 26);
			Ki_old_ = Ki[iout] = NV_Ith_S(dependent_variable__, 27);
			ato_f_old_ = ato_f[iout] = NV_Ith_S(dependent_variable__, 28);
			ito_f_old_ = ito_f[iout] = NV_Ith_S(dependent_variable__, 29);
			ato_s_old_ = ato_s[iout] = NV_Ith_S(dependent_variable__, 30);
			ito_s_old_ = ito_s[iout] = NV_Ith_S(dependent_variable__, 31);
			nKs_old_ = nKs[iout] = NV_Ith_S(dependent_variable__, 32);
			aur_old_ = aur[iout] = NV_Ith_S(dependent_variable__, 33);
			iur_old_ = iur[iout] = NV_Ith_S(dependent_variable__, 34);
			aKss_old_ = aKss[iout] = NV_Ith_S(dependent_variable__, 35);
			iKss_old_ = iKss[iout] = NV_Ith_S(dependent_variable__, 36);
			C_K2_old_ = C_K2[iout] = NV_Ith_S(dependent_variable__, 37);
			C_K1_old_ = C_K1[iout] = NV_Ith_S(dependent_variable__, 38);
			O_K_old_ = O_K[iout] = NV_Ith_S(dependent_variable__, 39);
			I_K_old_ = I_K[iout] = NV_Ith_S(dependent_variable__, 40);
                        istim[iout] = calc_Istim();

                        I_Ktof[iout] = calc_i_Kto_f();
                        I_Kur[iout] = calc_i_Kur();
                        I_Kss[iout] = calc_i_Kss();
                        I_Ks[iout] = calc_i_Ks();
                        I_Kr[iout] = calc_i_Kr();
                        I_K1[iout] = calc_i_K1();
                        I_CaL[iout] = calc_i_CaL();
                        I_Cab[iout] = calc_i_Cab();
                        I_Na[iout] = calc_i_Na();
                        I_NaK[iout] = calc_i_NaK();
                        I_NaCa[iout] = calc_i_NaCa();
                        I_Nab[iout] = calc_i_Nab();
                        I_pCa[iout] = calc_i_pCa();
                        I_ClCa[iout] = calc_i_ClCa();

			if (check_flag(&flag__, "CVode", 1)) break;
			if (flag__ == CV_SUCCESS){
				iout++;
				tout += dtime; // timestep
			}
			if (iout == num_iterations_bak ) break;
		}
	}
	double Solveode::solve(int firstcall__, int num_iterations__, int num_results__)
	{

		static int num_iterations_bak = 0;
		static int num_results_bak = 0;
		static int offset_step = 1;
		if(firstcall__){
			time_new = time;

			if(num_results__ <= 0)
				num_results__ = 1;
			if(num_iterations__ <= 0)
				num_iterations__ = 1;
			offset_step = num_iterations__ / num_results__;
			if(time_vec__ != NULL)free( time_vec__);
			time_vec__ = (double *)malloc(sizeof(double)*num_results__);
			V_old_ = V_ini_;
			if(V != NULL)free( V);
			V = (double *)malloc(sizeof(double)*num_results__);
			Cai_old_ = Cai_ini_;
			if(Cai != NULL)free( Cai);
			Cai = (double *)malloc(sizeof(double)*num_results__);
			Cass_old_ = Cass_ini_;
			if(Cass != NULL)free( Cass);
			Cass = (double *)malloc(sizeof(double)*num_results__);
			CaJSR_old_ = CaJSR_ini_;
			if(CaJSR != NULL)free( CaJSR);
			CaJSR = (double *)malloc(sizeof(double)*num_results__);
			CaNSR_old_ = CaNSR_ini_;
			if(CaNSR != NULL)free( CaNSR);
			CaNSR = (double *)malloc(sizeof(double)*num_results__);
			P_RyR_old_ = P_RyR_ini_;
			if(P_RyR != NULL)free( P_RyR);
			P_RyR = (double *)malloc(sizeof(double)*num_results__);
			LTRPN_Ca_old_ = LTRPN_Ca_ini_;
			if(LTRPN_Ca != NULL)free( LTRPN_Ca);
			LTRPN_Ca = (double *)malloc(sizeof(double)*num_results__);
			HTRPN_Ca_old_ = HTRPN_Ca_ini_;
			if(HTRPN_Ca != NULL)free( HTRPN_Ca);
			HTRPN_Ca = (double *)malloc(sizeof(double)*num_results__);
			P_O1_old_ = P_O1_ini_;
			if(P_O1 != NULL)free( P_O1);
			P_O1 = (double *)malloc(sizeof(double)*num_results__);
			P_O2_old_ = P_O2_ini_;
			if(P_O2 != NULL)free( P_O2);
			P_O2 = (double *)malloc(sizeof(double)*num_results__);
			P_C2_old_ = P_C2_ini_;
			if(P_C2 != NULL)free( P_C2);
			P_C2 = (double *)malloc(sizeof(double)*num_results__);
			O_old_ = O_ini_;
			if(O != NULL)free( O);
			O = (double *)malloc(sizeof(double)*num_results__);
			C2_old_ = C2_ini_;
			if(C2 != NULL)free( C2);
			C2 = (double *)malloc(sizeof(double)*num_results__);
			C3_old_ = C3_ini_;
			if(C3 != NULL)free( C3);
			C3 = (double *)malloc(sizeof(double)*num_results__);
			C4_old_ = C4_ini_;
			if(C4 != NULL)free( C4);
			C4 = (double *)malloc(sizeof(double)*num_results__);
			I1_old_ = I1_ini_;
			if(I1 != NULL)free( I1);
			I1 = (double *)malloc(sizeof(double)*num_results__);
			I2_old_ = I2_ini_;
			if(I2 != NULL)free( I2);
			I2 = (double *)malloc(sizeof(double)*num_results__);
			I3_old_ = I3_ini_;
			if(I3 != NULL)free( I3);
			I3 = (double *)malloc(sizeof(double)*num_results__);
			Nai_old_ = Nai_ini_;
			if(Nai != NULL)free( Nai);
			Nai = (double *)malloc(sizeof(double)*num_results__);
			C_Na2_old_ = C_Na2_ini_;
			if(C_Na2 != NULL)free( C_Na2);
			C_Na2 = (double *)malloc(sizeof(double)*num_results__);
			C_Na1_old_ = C_Na1_ini_;
			if(C_Na1 != NULL)free( C_Na1);
			C_Na1 = (double *)malloc(sizeof(double)*num_results__);
			O_Na_old_ = O_Na_ini_;
			if(O_Na != NULL)free( O_Na);
			O_Na = (double *)malloc(sizeof(double)*num_results__);
			IF_Na_old_ = IF_Na_ini_;
			if(IF_Na != NULL)free( IF_Na);
			IF_Na = (double *)malloc(sizeof(double)*num_results__);
			I1_Na_old_ = I1_Na_ini_;
			if(I1_Na != NULL)free( I1_Na);
			I1_Na = (double *)malloc(sizeof(double)*num_results__);
			I2_Na_old_ = I2_Na_ini_;
			if(I2_Na != NULL)free( I2_Na);
			I2_Na = (double *)malloc(sizeof(double)*num_results__);
			IC_Na2_old_ = IC_Na2_ini_;
			if(IC_Na2 != NULL)free( IC_Na2);
			IC_Na2 = (double *)malloc(sizeof(double)*num_results__);
			IC_Na3_old_ = IC_Na3_ini_;
			if(IC_Na3 != NULL)free( IC_Na3);
			IC_Na3 = (double *)malloc(sizeof(double)*num_results__);
			Ki_old_ = Ki_ini_;
			if(Ki != NULL)free( Ki);
			Ki = (double *)malloc(sizeof(double)*num_results__);
			ato_f_old_ = ato_f_ini_;
			if(ato_f != NULL)free( ato_f);
			ato_f = (double *)malloc(sizeof(double)*num_results__);
			ito_f_old_ = ito_f_ini_;
			if(ito_f != NULL)free( ito_f);
			ito_f = (double *)malloc(sizeof(double)*num_results__);
			ato_s_old_ = ato_s_ini_;
			if(ato_s != NULL)free( ato_s);
			ato_s = (double *)malloc(sizeof(double)*num_results__);
			ito_s_old_ = ito_s_ini_;
			if(ito_s != NULL)free( ito_s);
			ito_s = (double *)malloc(sizeof(double)*num_results__);
			nKs_old_ = nKs_ini_;
			if(nKs != NULL)free( nKs);
			nKs = (double *)malloc(sizeof(double)*num_results__);
			aur_old_ = aur_ini_;
			if(aur != NULL)free( aur);
			aur = (double *)malloc(sizeof(double)*num_results__);
			iur_old_ = iur_ini_;
			if(iur != NULL)free( iur);
			iur = (double *)malloc(sizeof(double)*num_results__);
			aKss_old_ = aKss_ini_;
			if(aKss != NULL)free( aKss);
			aKss = (double *)malloc(sizeof(double)*num_results__);
			iKss_old_ = iKss_ini_;
			if(iKss != NULL)free( iKss);
			iKss = (double *)malloc(sizeof(double)*num_results__);
			C_K2_old_ = C_K2_ini_;
			if(C_K2 != NULL)free( C_K2);
			C_K2 = (double *)malloc(sizeof(double)*num_results__);
			C_K1_old_ = C_K1_ini_;
			if(C_K1 != NULL)free( C_K1);
			C_K1 = (double *)malloc(sizeof(double)*num_results__);
			O_K_old_ = O_K_ini_;
			if(O_K != NULL)free( O_K);
			O_K = (double *)malloc(sizeof(double)*num_results__);
			I_K_old_ = I_K_ini_;
			if(I_K != NULL)free( I_K);
			I_K = (double *)malloc(sizeof(double)*num_results__);
			num_results_bak = num_results__;
			num_iterations_bak = num_iterations__;
		}
		int counter_it__ = 0;
		int aux = num_iterations_bak%num_results_bak;
		for(int it_countx = 1; it_countx<=num_iterations_bak; it_countx++){
			time_new += dtime;

			V_new_=dtime*(((-(calc_i_CaL()+calc_i_pCa()+calc_i_NaCa()+calc_i_Cab()+calc_i_Na()+calc_i_Nab()+calc_i_NaK()+calc_i_Kto_f()+calc_i_Kto_s()+calc_i_K1()+calc_i_Ks()+calc_i_Kur()+calc_i_Kss()+calc_i_Kr()+calc_i_ClCa()+calc_Istim()))/Cm))+V_old_;	// 0
			Cai_new_=dtime*((calc_Bi()*((calc_J_leak()+calc_J_xfer())-(calc_J_up()+calc_J_trpn()+((((calc_i_Cab()+calc_i_pCa())-(2*calc_i_NaCa()))*Acap*Cm)/(2*Vmyo*F))))))+Cai_old_;	// 1
			Cass_new_=dtime*((calc_Bss()*(((calc_J_rel()*VJSR)/Vss)-(((calc_J_xfer()*Vmyo)/Vss)+((calc_i_CaL()*Acap*Cm)/(2*Vss*F))))))+Cass_old_;	// 2
			CaJSR_new_=dtime*((calc_BJSR()*(calc_J_tr()-calc_J_rel())))+CaJSR_old_;	// 3
			CaNSR_new_=dtime*(((((calc_J_up()-calc_J_leak())*Vmyo)/VNSR)-((calc_J_tr()*VJSR)/VNSR)))+CaNSR_old_;	// 4
			P_RyR_new_=dtime*((((-0.04)*P_RyR_old_)-(((0.1*calc_i_CaL())/i_CaL_max)*exp(((-pow((V_old_-5),2))/648)))))+P_RyR_old_;	// 14
			LTRPN_Ca_new_=dtime*(((k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_))-(k_minus_ltrpn*LTRPN_Ca_old_)))+LTRPN_Ca_old_;	// 15
			HTRPN_Ca_new_=dtime*(((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))-(k_minus_htrpn*HTRPN_Ca_old_)))+HTRPN_Ca_old_;	// 16
			P_O1_new_=dtime*((((k_plus_a*pow(Cass_old_,n)*calc_P_C1())+(k_minus_b*P_O2_old_)+(k_minus_c*P_C2_old_))-((k_minus_a*P_O1_old_)+(k_plus_b*pow(Cass_old_,m)*P_O1_old_)+(k_plus_c*P_O1_old_))))+P_O1_old_;	// 17
			P_O2_new_=dtime*(((k_plus_b*pow(Cass_old_,m)*P_O1_old_)-(k_minus_b*P_O2_old_)))+P_O2_old_;	// 19
			P_C2_new_=dtime*(((k_plus_c*P_O1_old_)-(k_minus_c*P_C2_old_)))+P_C2_old_;	// 20
			O_new_=dtime*((((calc_alpha()*C4_old_)+(Kpcb*I1_old_)+(0.001*((calc_alpha()*I2_old_)-(calc_Kpcf()*O_old_))))-((4*calc_beta()*O_old_)+(calc_gamma()*O_old_))))+O_old_;	// 22
			C2_new_=dtime*((((4*calc_alpha()*calc_C1())+(2*calc_beta()*C3_old_))-((calc_beta()*C2_old_)+(3*calc_alpha()*C2_old_))))+C2_old_;	// 24
			C3_new_=dtime*((((3*calc_alpha()*C2_old_)+(3*calc_beta()*C4_old_))-((2*calc_beta()*C3_old_)+(2*calc_alpha()*C3_old_))))+C3_old_;	// 25
			C4_new_=dtime*((((2*calc_alpha()*C3_old_)+(4*calc_beta()*O_old_)+(0.01*((4*Kpcb*calc_beta()*I1_old_)-(calc_alpha()*calc_gamma()*C4_old_)))+(0.002*((4*calc_beta()*I2_old_)-(calc_Kpcf()*C4_old_)))+(4*calc_beta()*Kpcb*I3_old_))-((3*calc_beta()*C4_old_)+(calc_alpha()*C4_old_)+(calc_gamma()*calc_Kpcf()*C4_old_))))+C4_old_;	// 26
			I1_new_=dtime*((((calc_gamma()*O_old_)+(0.001*((calc_alpha()*I3_old_)-(calc_Kpcf()*I1_old_)))+(0.01*((calc_alpha()*calc_gamma()*C4_old_)-(4*calc_beta()*calc_Kpcf()*I1_old_))))-(Kpcb*I1_old_)))+I1_old_;	// 27
			I2_new_=dtime*((((0.001*((calc_Kpcf()*O_old_)-(calc_alpha()*I2_old_)))+(Kpcb*I3_old_)+(0.002*((calc_Kpcf()*C4_old_)-(4*calc_beta()*I2_old_))))-(calc_gamma()*I2_old_)))+I2_old_;	// 28
			I3_new_=dtime*((((0.001*((calc_Kpcf()*I1_old_)-(calc_alpha()*I3_old_)))+(calc_gamma()*I2_old_)+(calc_gamma()*calc_Kpcf()*C4_old_))-((4*calc_beta()*Kpcb*I3_old_)+(Kpcb*I3_old_))))+I3_old_;	// 29
			Nai_new_=dtime*((((-(calc_i_Na()+calc_i_Nab()+(3*calc_i_NaK())+(3*calc_i_NaCa())))*Acap*Cm)/(Vmyo*F)))+Nai_old_;	// 38
			C_Na2_new_=dtime*((((calc_alpha_Na11()*calc_C_Na3())+(calc_beta_Na12()*C_Na1_old_)+(calc_alpha_Na3()*IC_Na2_old_))-((calc_beta_Na11()*C_Na2_old_)+(calc_alpha_Na12()*C_Na2_old_)+(calc_beta_Na3()*C_Na2_old_))))+C_Na2_old_;	// 42
			C_Na1_new_=dtime*((((calc_alpha_Na12()*C_Na2_old_)+(calc_beta_Na13()*O_Na_old_)+(calc_alpha_Na3()*IF_Na_old_))-((calc_beta_Na12()*C_Na1_old_)+(calc_alpha_Na13()*C_Na1_old_)+(calc_beta_Na3()*C_Na1_old_))))+C_Na1_old_;	// 43
			O_Na_new_=dtime*((((calc_alpha_Na13()*C_Na1_old_)+(calc_beta_Na2()*IF_Na_old_))-((calc_beta_Na13()*O_Na_old_)+(calc_alpha_Na2()*O_Na_old_))))+O_Na_old_;	// 44
			IF_Na_new_=dtime*((((calc_alpha_Na2()*O_Na_old_)+(calc_beta_Na3()*C_Na1_old_)+(calc_beta_Na4()*I1_Na_old_)+(calc_alpha_Na12()*IC_Na2_old_))-((calc_beta_Na2()*IF_Na_old_)+(calc_alpha_Na3()*IF_Na_old_)+(calc_alpha_Na4()*IF_Na_old_)+(calc_beta_Na12()*IF_Na_old_))))+IF_Na_old_;	// 45
			I1_Na_new_=dtime*((((calc_alpha_Na4()*IF_Na_old_)+(calc_beta_Na5()*I2_Na_old_))-((calc_beta_Na4()*I1_Na_old_)+(calc_alpha_Na5()*I1_Na_old_))))+I1_Na_old_;	// 46
			I2_Na_new_=dtime*(((calc_alpha_Na5()*I1_Na_old_)-(calc_beta_Na5()*I2_Na_old_)))+I2_Na_old_;	// 47
			IC_Na2_new_=dtime*((((calc_alpha_Na11()*IC_Na3_old_)+(calc_beta_Na12()*IF_Na_old_)+(calc_beta_Na3()*IC_Na2_old_))-((calc_beta_Na11()*IC_Na2_old_)+(calc_alpha_Na12()*IC_Na2_old_)+(calc_alpha_Na3()*IC_Na2_old_))))+IC_Na2_old_;	// 48
			IC_Na3_new_=dtime*((((calc_beta_Na11()*IC_Na2_old_)+(calc_beta_Na3()*calc_C_Na3()))-((calc_alpha_Na11()*IC_Na3_old_)+(calc_alpha_Na3()*IC_Na3_old_))))+IC_Na3_old_;	// 49
			Ki_new_=dtime*((((-((calc_i_Kto_f()+calc_i_Kto_s()+calc_i_K1()+calc_i_Ks()+calc_i_Kss()+calc_i_Kur()+calc_i_Kr())-(2*calc_i_NaK())))*Acap*Cm)/(Vmyo*F)))+Ki_old_;	// 65
			ato_f_new_=dtime*(((calc_alpha_a()*(1-ato_f_old_))-(calc_beta_a()*ato_f_old_)))+ato_f_old_;	// 68
			ito_f_new_=dtime*(((calc_alpha_i()*(1-ito_f_old_))-(calc_beta_i()*ito_f_old_)))+ito_f_old_;	// 69
			ato_s_new_=dtime*(((calc_ass()-ato_s_old_)/calc_tau_ta_s()))+ato_s_old_;	// 75
			ito_s_new_=dtime*(((calc_iss()-ito_s_old_)/calc_tau_ti_s()))+ito_s_old_;	// 76
			nKs_new_=dtime*(((calc_alpha_n()*(1-nKs_old_))-(calc_beta_n()*nKs_old_)))+nKs_old_;	// 83
			aur_new_=dtime*(((calc_ass()-aur_old_)/calc_tau_aur()))+aur_old_;	// 87
			iur_new_=dtime*(((calc_iss()-iur_old_)/calc_tau_iur()))+iur_old_;	// 88
			aKss_new_=dtime*(((calc_ass()-aKss_old_)/calc_tau_Kss()))+aKss_old_;	// 92
			iKss_new_=dtime*(0)+iKss_old_;	// 93
			C_K2_new_=dtime*((((kf*C_K1_old_)+(calc_beta_a1()*O_K_old_))-((kb*C_K2_old_)+(calc_alpha_a1()*C_K2_old_))))+C_K2_old_;	// 97
			C_K1_new_=dtime*((((calc_alpha_a0()*calc_C_K0())+(kb*C_K2_old_))-((calc_beta_a0()*C_K1_old_)+(kf*C_K1_old_))))+C_K1_old_;	// 98
			O_K_new_=dtime*((((calc_alpha_a1()*C_K2_old_)+(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current()*I_K_old_))-((calc_beta_a1()*O_K_old_)+(calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current()*O_K_old_))))+O_K_old_;	// 99
			I_K_new_=dtime*(((calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current()*O_K_old_)-(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current()*I_K_old_)))+I_K_old_;	// 100

			if(it_countx != aux && (it_countx-aux)%offset_step == 0){
				V[counter_it__] = V_new_;
				Cai[counter_it__] = Cai_new_;
				Cass[counter_it__] = Cass_new_;
				CaJSR[counter_it__] = CaJSR_new_;
				CaNSR[counter_it__] = CaNSR_new_;
				P_RyR[counter_it__] = P_RyR_new_;
				LTRPN_Ca[counter_it__] = LTRPN_Ca_new_;
				HTRPN_Ca[counter_it__] = HTRPN_Ca_new_;
				P_O1[counter_it__] = P_O1_new_;
				P_O2[counter_it__] = P_O2_new_;
				P_C2[counter_it__] = P_C2_new_;
				O[counter_it__] = O_new_;
				C2[counter_it__] = C2_new_;
				C3[counter_it__] = C3_new_;
				C4[counter_it__] = C4_new_;
				I1[counter_it__] = I1_new_;
				I2[counter_it__] = I2_new_;
				I3[counter_it__] = I3_new_;
				Nai[counter_it__] = Nai_new_;
				C_Na2[counter_it__] = C_Na2_new_;
				C_Na1[counter_it__] = C_Na1_new_;
				O_Na[counter_it__] = O_Na_new_;
				IF_Na[counter_it__] = IF_Na_new_;
				I1_Na[counter_it__] = I1_Na_new_;
				I2_Na[counter_it__] = I2_Na_new_;
				IC_Na2[counter_it__] = IC_Na2_new_;
				IC_Na3[counter_it__] = IC_Na3_new_;
				Ki[counter_it__] = Ki_new_;
				ato_f[counter_it__] = ato_f_new_;
				ito_f[counter_it__] = ito_f_new_;
				ato_s[counter_it__] = ato_s_new_;
				ito_s[counter_it__] = ito_s_new_;
				nKs[counter_it__] = nKs_new_;
				aur[counter_it__] = aur_new_;
				iur[counter_it__] = iur_new_;
				aKss[counter_it__] = aKss_new_;
				iKss[counter_it__] = iKss_new_;
				C_K2[counter_it__] = C_K2_new_;
				C_K1[counter_it__] = C_K1_new_;
				O_K[counter_it__] = O_K_new_;
				I_K[counter_it__] = I_K_new_;
				time_vec__[counter_it__] = time_new;
				counter_it__++;
			}
			V_old_ = V_new_;
			Cai_old_ = Cai_new_;
			Cass_old_ = Cass_new_;
			CaJSR_old_ = CaJSR_new_;
			CaNSR_old_ = CaNSR_new_;
			P_RyR_old_ = P_RyR_new_;
			LTRPN_Ca_old_ = LTRPN_Ca_new_;
			HTRPN_Ca_old_ = HTRPN_Ca_new_;
			P_O1_old_ = P_O1_new_;
			P_O2_old_ = P_O2_new_;
			P_C2_old_ = P_C2_new_;
			O_old_ = O_new_;
			C2_old_ = C2_new_;
			C3_old_ = C3_new_;
			C4_old_ = C4_new_;
			I1_old_ = I1_new_;
			I2_old_ = I2_new_;
			I3_old_ = I3_new_;
			Nai_old_ = Nai_new_;
			C_Na2_old_ = C_Na2_new_;
			C_Na1_old_ = C_Na1_new_;
			O_Na_old_ = O_Na_new_;
			IF_Na_old_ = IF_Na_new_;
			I1_Na_old_ = I1_Na_new_;
			I2_Na_old_ = I2_Na_new_;
			IC_Na2_old_ = IC_Na2_new_;
			IC_Na3_old_ = IC_Na3_new_;
			Ki_old_ = Ki_new_;
			ato_f_old_ = ato_f_new_;
			ito_f_old_ = ito_f_new_;
			ato_s_old_ = ato_s_new_;
			ito_s_old_ = ito_s_new_;
			nKs_old_ = nKs_new_;
			aur_old_ = aur_new_;
			iur_old_ = iur_new_;
			aKss_old_ = aKss_new_;
			iKss_old_ = iKss_new_;
			C_K2_old_ = C_K2_new_;
			C_K1_old_ = C_K1_new_;
			O_K_old_ = O_K_new_;
			I_K_old_ = I_K_new_;
		}
		return (num_iterations_bak%offset_step)*dtime;
	}

	double* Solveode::getIndependentVar()
	{
		return time_vec__;
	}

	double* Solveode::getSolution(int indVariable)
	{
		switch(indVariable)
		{
		case 0:		return V;    break;
		case 1:		return Cai;    break;
		case 2:		return Cass;    break;
		case 3:		return CaJSR;    break;
		case 4:		return CaNSR;    break;
		case 5:		return P_RyR;    break;
		case 6:		return LTRPN_Ca;    break;
		case 7:		return HTRPN_Ca;    break;
		case 8:		return P_O1;    break;
		case 9:		return P_O2;    break;
		case 10:		return P_C2;    break;
		case 11:		return O;    break;
		case 12:		return C2;    break;
		case 13:		return C3;    break;
		case 14:		return C4;    break;
		case 15:		return I1;    break;
		case 16:		return I2;    break;
		case 17:		return I3;    break;
		case 18:		return Nai;    break;
		case 19:		return C_Na2;    break;
		case 20:		return C_Na1;    break;
		case 21:		return O_Na;    break;
		case 22:		return IF_Na;    break;
		case 23:		return I1_Na;    break;
		case 24:		return I2_Na;    break;
		case 25:		return IC_Na2;    break;
		case 26:		return IC_Na3;    break;
		case 27:		return Ki;    break;
		case 28:		return ato_f;    break;
		case 29:		return ito_f;    break;
		case 30:		return ato_s;    break;
		case 31:		return ito_s;    break;
		case 32:		return nKs;    break;
		case 33:		return aur;    break;
		case 34:		return iur;    break;
		case 35:		return aKss;    break;
		case 36:		return iKss;    break;
		case 37:		return C_K2;    break;
		case 38:		return C_K1;    break;
		case 39:		return O_K;    break;
		case 40:		return I_K;    break;
                case 41:                return istim;    break;
                case 100:                return I_Ktof;    break;
                case 101:                return I_Kur;    break;
                case 102:                return I_Kss;    break;
                case 103:                return I_Ks;    break;
                case 104:                return I_Kr;    break;
                case 105:                return I_K1;    break;
                case 106:                return I_CaL;    break;
                case 107:                return I_Cab;    break;
                case 108:                return I_Na;    break;
                case 109:                return I_NaK;    break;
                case 110:                return I_NaCa;    break;
                case 111:                return I_Nab;    break;
                case 112:                return I_pCa;    break;
                case 113:                return I_ClCa;    break;
		default:	return NULL;    break;
		}
	}
	double Solveode::solveToFile(char *filename, char *fileaccess, int firstcall__, int num_iterations__, int num_results__)
	{

		static int num_iterations_bak = 0;
		static int num_results_bak = 0;
		static int offset_step = 1;
		static char *fileaccess_bak = "";
		if(firstcall__){
			time_new = time;

			if(num_results__ <= 0)
				num_results__ = 1;
			if(num_iterations__ <= 0)
				num_iterations__ = 1;
			offset_step = num_iterations__ / num_results__;
			V_old_ = V_ini_;
			if(V != NULL)free( V);
			V = (double *)malloc(sizeof(double)*num_results__);
			Cai_old_ = Cai_ini_;
			if(Cai != NULL)free( Cai);
			Cai = (double *)malloc(sizeof(double)*num_results__);
			Cass_old_ = Cass_ini_;
			if(Cass != NULL)free( Cass);
			Cass = (double *)malloc(sizeof(double)*num_results__);
			CaJSR_old_ = CaJSR_ini_;
			if(CaJSR != NULL)free( CaJSR);
			CaJSR = (double *)malloc(sizeof(double)*num_results__);
			CaNSR_old_ = CaNSR_ini_;
			if(CaNSR != NULL)free( CaNSR);
			CaNSR = (double *)malloc(sizeof(double)*num_results__);
			P_RyR_old_ = P_RyR_ini_;
			if(P_RyR != NULL)free( P_RyR);
			P_RyR = (double *)malloc(sizeof(double)*num_results__);
			LTRPN_Ca_old_ = LTRPN_Ca_ini_;
			if(LTRPN_Ca != NULL)free( LTRPN_Ca);
			LTRPN_Ca = (double *)malloc(sizeof(double)*num_results__);
			HTRPN_Ca_old_ = HTRPN_Ca_ini_;
			if(HTRPN_Ca != NULL)free( HTRPN_Ca);
			HTRPN_Ca = (double *)malloc(sizeof(double)*num_results__);
			P_O1_old_ = P_O1_ini_;
			if(P_O1 != NULL)free( P_O1);
			P_O1 = (double *)malloc(sizeof(double)*num_results__);
			P_O2_old_ = P_O2_ini_;
			if(P_O2 != NULL)free( P_O2);
			P_O2 = (double *)malloc(sizeof(double)*num_results__);
			P_C2_old_ = P_C2_ini_;
			if(P_C2 != NULL)free( P_C2);
			P_C2 = (double *)malloc(sizeof(double)*num_results__);
			O_old_ = O_ini_;
			if(O != NULL)free( O);
			O = (double *)malloc(sizeof(double)*num_results__);
			C2_old_ = C2_ini_;
			if(C2 != NULL)free( C2);
			C2 = (double *)malloc(sizeof(double)*num_results__);
			C3_old_ = C3_ini_;
			if(C3 != NULL)free( C3);
			C3 = (double *)malloc(sizeof(double)*num_results__);
			C4_old_ = C4_ini_;
			if(C4 != NULL)free( C4);
			C4 = (double *)malloc(sizeof(double)*num_results__);
			I1_old_ = I1_ini_;
			if(I1 != NULL)free( I1);
			I1 = (double *)malloc(sizeof(double)*num_results__);
			I2_old_ = I2_ini_;
			if(I2 != NULL)free( I2);
			I2 = (double *)malloc(sizeof(double)*num_results__);
			I3_old_ = I3_ini_;
			if(I3 != NULL)free( I3);
			I3 = (double *)malloc(sizeof(double)*num_results__);
			Nai_old_ = Nai_ini_;
			if(Nai != NULL)free( Nai);
			Nai = (double *)malloc(sizeof(double)*num_results__);
			C_Na2_old_ = C_Na2_ini_;
			if(C_Na2 != NULL)free( C_Na2);
			C_Na2 = (double *)malloc(sizeof(double)*num_results__);
			C_Na1_old_ = C_Na1_ini_;
			if(C_Na1 != NULL)free( C_Na1);
			C_Na1 = (double *)malloc(sizeof(double)*num_results__);
			O_Na_old_ = O_Na_ini_;
			if(O_Na != NULL)free( O_Na);
			O_Na = (double *)malloc(sizeof(double)*num_results__);
			IF_Na_old_ = IF_Na_ini_;
			if(IF_Na != NULL)free( IF_Na);
			IF_Na = (double *)malloc(sizeof(double)*num_results__);
			I1_Na_old_ = I1_Na_ini_;
			if(I1_Na != NULL)free( I1_Na);
			I1_Na = (double *)malloc(sizeof(double)*num_results__);
			I2_Na_old_ = I2_Na_ini_;
			if(I2_Na != NULL)free( I2_Na);
			I2_Na = (double *)malloc(sizeof(double)*num_results__);
			IC_Na2_old_ = IC_Na2_ini_;
			if(IC_Na2 != NULL)free( IC_Na2);
			IC_Na2 = (double *)malloc(sizeof(double)*num_results__);
			IC_Na3_old_ = IC_Na3_ini_;
			if(IC_Na3 != NULL)free( IC_Na3);
			IC_Na3 = (double *)malloc(sizeof(double)*num_results__);
			Ki_old_ = Ki_ini_;
			if(Ki != NULL)free( Ki);
			Ki = (double *)malloc(sizeof(double)*num_results__);
			ato_f_old_ = ato_f_ini_;
			if(ato_f != NULL)free( ato_f);
			ato_f = (double *)malloc(sizeof(double)*num_results__);
			ito_f_old_ = ito_f_ini_;
			if(ito_f != NULL)free( ito_f);
			ito_f = (double *)malloc(sizeof(double)*num_results__);
			ato_s_old_ = ato_s_ini_;
			if(ato_s != NULL)free( ato_s);
			ato_s = (double *)malloc(sizeof(double)*num_results__);
			ito_s_old_ = ito_s_ini_;
			if(ito_s != NULL)free( ito_s);
			ito_s = (double *)malloc(sizeof(double)*num_results__);
			nKs_old_ = nKs_ini_;
			if(nKs != NULL)free( nKs);
			nKs = (double *)malloc(sizeof(double)*num_results__);
			aur_old_ = aur_ini_;
			if(aur != NULL)free( aur);
			aur = (double *)malloc(sizeof(double)*num_results__);
			iur_old_ = iur_ini_;
			if(iur != NULL)free( iur);
			iur = (double *)malloc(sizeof(double)*num_results__);
			aKss_old_ = aKss_ini_;
			if(aKss != NULL)free( aKss);
			aKss = (double *)malloc(sizeof(double)*num_results__);
			iKss_old_ = iKss_ini_;
			if(iKss != NULL)free( iKss);
			iKss = (double *)malloc(sizeof(double)*num_results__);
			C_K2_old_ = C_K2_ini_;
			if(C_K2 != NULL)free( C_K2);
			C_K2 = (double *)malloc(sizeof(double)*num_results__);
			C_K1_old_ = C_K1_ini_;
			if(C_K1 != NULL)free( C_K1);
			C_K1 = (double *)malloc(sizeof(double)*num_results__);
			O_K_old_ = O_K_ini_;
			if(O_K != NULL)free( O_K);
			O_K = (double *)malloc(sizeof(double)*num_results__);
			I_K_old_ = I_K_ini_;
			if(I_K != NULL)free( I_K);
			I_K = (double *)malloc(sizeof(double)*num_results__);
			num_results_bak = num_results__;
			num_iterations_bak = num_iterations__;
			fileaccess_bak = fileaccess;
		}
		FILE *file = fopen(filename, fileaccess_bak);
		if(!file){
			fprintf(stderr,"ERROR - solveToFile - Unable to open file %s\n",filename);
			exit(1);
		}
		int counter_it__ = 0;
		int aux = num_iterations_bak%num_results_bak;
		for(int it_countx = 1; it_countx<=num_iterations_bak; it_countx++){
			time_new += dtime;

			V_new_=dtime*(((-(calc_i_CaL()+calc_i_pCa()+calc_i_NaCa()+calc_i_Cab()+calc_i_Na()+calc_i_Nab()+calc_i_NaK()+calc_i_Kto_f()+calc_i_Kto_s()+calc_i_K1()+calc_i_Ks()+calc_i_Kur()+calc_i_Kss()+calc_i_Kr()+calc_i_ClCa()+calc_Istim()))/Cm))+V_old_;	// 0
			Cai_new_=dtime*((calc_Bi()*((calc_J_leak()+calc_J_xfer())-(calc_J_up()+calc_J_trpn()+((((calc_i_Cab()+calc_i_pCa())-(2*calc_i_NaCa()))*Acap*Cm)/(2*Vmyo*F))))))+Cai_old_;	// 1
			Cass_new_=dtime*((calc_Bss()*(((calc_J_rel()*VJSR)/Vss)-(((calc_J_xfer()*Vmyo)/Vss)+((calc_i_CaL()*Acap*Cm)/(2*Vss*F))))))+Cass_old_;	// 2
			CaJSR_new_=dtime*((calc_BJSR()*(calc_J_tr()-calc_J_rel())))+CaJSR_old_;	// 3
			CaNSR_new_=dtime*(((((calc_J_up()-calc_J_leak())*Vmyo)/VNSR)-((calc_J_tr()*VJSR)/VNSR)))+CaNSR_old_;	// 4
			P_RyR_new_=dtime*((((-0.04)*P_RyR_old_)-(((0.1*calc_i_CaL())/i_CaL_max)*exp(((-pow((V_old_-5),2))/648)))))+P_RyR_old_;	// 14
			LTRPN_Ca_new_=dtime*(((k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_))-(k_minus_ltrpn*LTRPN_Ca_old_)))+LTRPN_Ca_old_;	// 15
			HTRPN_Ca_new_=dtime*(((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))-(k_minus_htrpn*HTRPN_Ca_old_)))+HTRPN_Ca_old_;	// 16
			P_O1_new_=dtime*((((k_plus_a*pow(Cass_old_,n)*calc_P_C1())+(k_minus_b*P_O2_old_)+(k_minus_c*P_C2_old_))-((k_minus_a*P_O1_old_)+(k_plus_b*pow(Cass_old_,m)*P_O1_old_)+(k_plus_c*P_O1_old_))))+P_O1_old_;	// 17
			P_O2_new_=dtime*(((k_plus_b*pow(Cass_old_,m)*P_O1_old_)-(k_minus_b*P_O2_old_)))+P_O2_old_;	// 19
			P_C2_new_=dtime*(((k_plus_c*P_O1_old_)-(k_minus_c*P_C2_old_)))+P_C2_old_;	// 20
			O_new_=dtime*((((calc_alpha()*C4_old_)+(Kpcb*I1_old_)+(0.001*((calc_alpha()*I2_old_)-(calc_Kpcf()*O_old_))))-((4*calc_beta()*O_old_)+(calc_gamma()*O_old_))))+O_old_;	// 22
			C2_new_=dtime*((((4*calc_alpha()*calc_C1())+(2*calc_beta()*C3_old_))-((calc_beta()*C2_old_)+(3*calc_alpha()*C2_old_))))+C2_old_;	// 24
			C3_new_=dtime*((((3*calc_alpha()*C2_old_)+(3*calc_beta()*C4_old_))-((2*calc_beta()*C3_old_)+(2*calc_alpha()*C3_old_))))+C3_old_;	// 25
			C4_new_=dtime*((((2*calc_alpha()*C3_old_)+(4*calc_beta()*O_old_)+(0.01*((4*Kpcb*calc_beta()*I1_old_)-(calc_alpha()*calc_gamma()*C4_old_)))+(0.002*((4*calc_beta()*I2_old_)-(calc_Kpcf()*C4_old_)))+(4*calc_beta()*Kpcb*I3_old_))-((3*calc_beta()*C4_old_)+(calc_alpha()*C4_old_)+(calc_gamma()*calc_Kpcf()*C4_old_))))+C4_old_;	// 26
			I1_new_=dtime*((((calc_gamma()*O_old_)+(0.001*((calc_alpha()*I3_old_)-(calc_Kpcf()*I1_old_)))+(0.01*((calc_alpha()*calc_gamma()*C4_old_)-(4*calc_beta()*calc_Kpcf()*I1_old_))))-(Kpcb*I1_old_)))+I1_old_;	// 27
			I2_new_=dtime*((((0.001*((calc_Kpcf()*O_old_)-(calc_alpha()*I2_old_)))+(Kpcb*I3_old_)+(0.002*((calc_Kpcf()*C4_old_)-(4*calc_beta()*I2_old_))))-(calc_gamma()*I2_old_)))+I2_old_;	// 28
			I3_new_=dtime*((((0.001*((calc_Kpcf()*I1_old_)-(calc_alpha()*I3_old_)))+(calc_gamma()*I2_old_)+(calc_gamma()*calc_Kpcf()*C4_old_))-((4*calc_beta()*Kpcb*I3_old_)+(Kpcb*I3_old_))))+I3_old_;	// 29
			Nai_new_=dtime*((((-(calc_i_Na()+calc_i_Nab()+(3*calc_i_NaK())+(3*calc_i_NaCa())))*Acap*Cm)/(Vmyo*F)))+Nai_old_;	// 38
			C_Na2_new_=dtime*((((calc_alpha_Na11()*calc_C_Na3())+(calc_beta_Na12()*C_Na1_old_)+(calc_alpha_Na3()*IC_Na2_old_))-((calc_beta_Na11()*C_Na2_old_)+(calc_alpha_Na12()*C_Na2_old_)+(calc_beta_Na3()*C_Na2_old_))))+C_Na2_old_;	// 42
			C_Na1_new_=dtime*((((calc_alpha_Na12()*C_Na2_old_)+(calc_beta_Na13()*O_Na_old_)+(calc_alpha_Na3()*IF_Na_old_))-((calc_beta_Na12()*C_Na1_old_)+(calc_alpha_Na13()*C_Na1_old_)+(calc_beta_Na3()*C_Na1_old_))))+C_Na1_old_;	// 43
			O_Na_new_=dtime*((((calc_alpha_Na13()*C_Na1_old_)+(calc_beta_Na2()*IF_Na_old_))-((calc_beta_Na13()*O_Na_old_)+(calc_alpha_Na2()*O_Na_old_))))+O_Na_old_;	// 44
			IF_Na_new_=dtime*((((calc_alpha_Na2()*O_Na_old_)+(calc_beta_Na3()*C_Na1_old_)+(calc_beta_Na4()*I1_Na_old_)+(calc_alpha_Na12()*IC_Na2_old_))-((calc_beta_Na2()*IF_Na_old_)+(calc_alpha_Na3()*IF_Na_old_)+(calc_alpha_Na4()*IF_Na_old_)+(calc_beta_Na12()*IF_Na_old_))))+IF_Na_old_;	// 45
			I1_Na_new_=dtime*((((calc_alpha_Na4()*IF_Na_old_)+(calc_beta_Na5()*I2_Na_old_))-((calc_beta_Na4()*I1_Na_old_)+(calc_alpha_Na5()*I1_Na_old_))))+I1_Na_old_;	// 46
			I2_Na_new_=dtime*(((calc_alpha_Na5()*I1_Na_old_)-(calc_beta_Na5()*I2_Na_old_)))+I2_Na_old_;	// 47
			IC_Na2_new_=dtime*((((calc_alpha_Na11()*IC_Na3_old_)+(calc_beta_Na12()*IF_Na_old_)+(calc_beta_Na3()*IC_Na2_old_))-((calc_beta_Na11()*IC_Na2_old_)+(calc_alpha_Na12()*IC_Na2_old_)+(calc_alpha_Na3()*IC_Na2_old_))))+IC_Na2_old_;	// 48
			IC_Na3_new_=dtime*((((calc_beta_Na11()*IC_Na2_old_)+(calc_beta_Na3()*calc_C_Na3()))-((calc_alpha_Na11()*IC_Na3_old_)+(calc_alpha_Na3()*IC_Na3_old_))))+IC_Na3_old_;	// 49
			Ki_new_=dtime*((((-((calc_i_Kto_f()+calc_i_Kto_s()+calc_i_K1()+calc_i_Ks()+calc_i_Kss()+calc_i_Kur()+calc_i_Kr())-(2*calc_i_NaK())))*Acap*Cm)/(Vmyo*F)))+Ki_old_;	// 65
			ato_f_new_=dtime*(((calc_alpha_a()*(1-ato_f_old_))-(calc_beta_a()*ato_f_old_)))+ato_f_old_;	// 68
			ito_f_new_=dtime*(((calc_alpha_i()*(1-ito_f_old_))-(calc_beta_i()*ito_f_old_)))+ito_f_old_;	// 69
			ato_s_new_=dtime*(((calc_ass()-ato_s_old_)/calc_tau_ta_s()))+ato_s_old_;	// 75
			ito_s_new_=dtime*(((calc_iss()-ito_s_old_)/calc_tau_ti_s()))+ito_s_old_;	// 76
			nKs_new_=dtime*(((calc_alpha_n()*(1-nKs_old_))-(calc_beta_n()*nKs_old_)))+nKs_old_;	// 83
			aur_new_=dtime*(((calc_ass()-aur_old_)/calc_tau_aur()))+aur_old_;	// 87
			iur_new_=dtime*(((calc_iss()-iur_old_)/calc_tau_iur()))+iur_old_;	// 88
			aKss_new_=dtime*(((calc_ass()-aKss_old_)/calc_tau_Kss()))+aKss_old_;	// 92
			iKss_new_=dtime*(0)+iKss_old_;	// 93
			C_K2_new_=dtime*((((kf*C_K1_old_)+(calc_beta_a1()*O_K_old_))-((kb*C_K2_old_)+(calc_alpha_a1()*C_K2_old_))))+C_K2_old_;	// 97
			C_K1_new_=dtime*((((calc_alpha_a0()*calc_C_K0())+(kb*C_K2_old_))-((calc_beta_a0()*C_K1_old_)+(kf*C_K1_old_))))+C_K1_old_;	// 98
			O_K_new_=dtime*((((calc_alpha_a1()*C_K2_old_)+(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current()*I_K_old_))-((calc_beta_a1()*O_K_old_)+(calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current()*O_K_old_))))+O_K_old_;	// 99
			I_K_new_=dtime*(((calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current()*O_K_old_)-(calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current()*I_K_old_)))+I_K_old_;	// 100

			if(it_countx != aux && (it_countx-aux)%offset_step == 0){
				fprintf(file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",time_new,V_new_,Cai_new_,Cass_new_,CaJSR_new_,CaNSR_new_,P_RyR_new_,LTRPN_Ca_new_,HTRPN_Ca_new_,P_O1_new_,P_O2_new_,P_C2_new_,O_new_,C2_new_,C3_new_,C4_new_,I1_new_,I2_new_,I3_new_,Nai_new_,C_Na2_new_,C_Na1_new_,O_Na_new_,IF_Na_new_,I1_Na_new_,I2_Na_new_,IC_Na2_new_,IC_Na3_new_,Ki_new_,ato_f_new_,ito_f_new_,ato_s_new_,ito_s_new_,nKs_new_,aur_new_,iur_new_,aKss_new_,iKss_new_,C_K2_new_,C_K1_new_,O_K_new_,I_K_new_);
				counter_it__++;
			}
			V_old_ = V_new_;
			Cai_old_ = Cai_new_;
			Cass_old_ = Cass_new_;
			CaJSR_old_ = CaJSR_new_;
			CaNSR_old_ = CaNSR_new_;
			P_RyR_old_ = P_RyR_new_;
			LTRPN_Ca_old_ = LTRPN_Ca_new_;
			HTRPN_Ca_old_ = HTRPN_Ca_new_;
			P_O1_old_ = P_O1_new_;
			P_O2_old_ = P_O2_new_;
			P_C2_old_ = P_C2_new_;
			O_old_ = O_new_;
			C2_old_ = C2_new_;
			C3_old_ = C3_new_;
			C4_old_ = C4_new_;
			I1_old_ = I1_new_;
			I2_old_ = I2_new_;
			I3_old_ = I3_new_;
			Nai_old_ = Nai_new_;
			C_Na2_old_ = C_Na2_new_;
			C_Na1_old_ = C_Na1_new_;
			O_Na_old_ = O_Na_new_;
			IF_Na_old_ = IF_Na_new_;
			I1_Na_old_ = I1_Na_new_;
			I2_Na_old_ = I2_Na_new_;
			IC_Na2_old_ = IC_Na2_new_;
			IC_Na3_old_ = IC_Na3_new_;
			Ki_old_ = Ki_new_;
			ato_f_old_ = ato_f_new_;
			ito_f_old_ = ito_f_new_;
			ato_s_old_ = ato_s_new_;
			ito_s_old_ = ito_s_new_;
			nKs_old_ = nKs_new_;
			aur_old_ = aur_new_;
			iur_old_ = iur_new_;
			aKss_old_ = aKss_new_;
			iKss_old_ = iKss_new_;
			C_K2_old_ = C_K2_new_;
			C_K1_old_ = C_K1_new_;
			O_K_old_ = O_K_new_;
			I_K_old_ = I_K_new_;
		}
		fclose(file);
		return (num_iterations_bak%offset_step)*dtime;
	}
	void Solveode::solveCVODEToFile(char *filename, char *fileaccess, int firstcall__, int steps__)
	{

		static int num_iterations_bak = 0;
		static char *fileaccess_bak = "";
		if(firstcall__){
			if(steps__ <= 0)
				steps__ = 1;
			num_iterations_bak = steps__;
			fileaccess_bak = fileaccess;
		}
		FILE *file = fopen(filename, fileaccess_bak);
		if(!file){
			fprintf(stderr,"ERROR - solveCVODEToFile - Unable to open file %s\n",filename);
			exit(1);
		}
		solveCVODE(firstcall__, num_iterations_bak);
		for(int i=0;i<num_iterations_bak; i++){
			fprintf(file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",time_vec__[i],V[i],Cai[i],Cass[i],CaJSR[i],CaNSR[i],P_RyR[i],LTRPN_Ca[i],HTRPN_Ca[i],P_O1[i],P_O2[i],P_C2[i],O[i],C2[i],C3[i],C4[i],I1[i],I2[i],I3[i],Nai[i],C_Na2[i],C_Na1[i],O_Na[i],IF_Na[i],I1_Na[i],I2_Na[i],IC_Na2[i],IC_Na3[i],Ki[i],ato_f[i],ito_f[i],ato_s[i],ito_s[i],nKs[i],aur[i],iur[i],aKss[i],iKss[i],C_K2[i],C_K1[i],O_K[i],I_K[i]);
		}
		fclose(file);
	}

	void Solveode::reInitCVODE()
	{

		flag__ = CVodeReInit(cvode_mem_cvode__, f__, time, dependent_variable__, CV_SV, reltol__, abstol__);
		if (check_flag(&flag__, "CVodeReInit", 1))
			exit(1);
	}

	void Solveode::setReltol(double value_new)
	{
		reltol__ = value_new;
	}

	void Solveode::setAbstol(int index, double value_new)
	{
		switch(index){
		case 0:		NV_Ith_S(abstol__, 0) = value_new;    break;
		case 1:		NV_Ith_S(abstol__, 1) = value_new;    break;
		case 2:		NV_Ith_S(abstol__, 2) = value_new;    break;
		case 3:		NV_Ith_S(abstol__, 3) = value_new;    break;
		case 4:		NV_Ith_S(abstol__, 4) = value_new;    break;
		case 5:		NV_Ith_S(abstol__, 5) = value_new;    break;
		case 6:		NV_Ith_S(abstol__, 6) = value_new;    break;
		case 7:		NV_Ith_S(abstol__, 7) = value_new;    break;
		case 8:		NV_Ith_S(abstol__, 8) = value_new;    break;
		case 9:		NV_Ith_S(abstol__, 9) = value_new;    break;
		case 10:		NV_Ith_S(abstol__, 10) = value_new;    break;
		case 11:		NV_Ith_S(abstol__, 11) = value_new;    break;
		case 12:		NV_Ith_S(abstol__, 12) = value_new;    break;
		case 13:		NV_Ith_S(abstol__, 13) = value_new;    break;
		case 14:		NV_Ith_S(abstol__, 14) = value_new;    break;
		case 15:		NV_Ith_S(abstol__, 15) = value_new;    break;
		case 16:		NV_Ith_S(abstol__, 16) = value_new;    break;
		case 17:		NV_Ith_S(abstol__, 17) = value_new;    break;
		case 18:		NV_Ith_S(abstol__, 18) = value_new;    break;
		case 19:		NV_Ith_S(abstol__, 19) = value_new;    break;
		case 20:		NV_Ith_S(abstol__, 20) = value_new;    break;
		case 21:		NV_Ith_S(abstol__, 21) = value_new;    break;
		case 22:		NV_Ith_S(abstol__, 22) = value_new;    break;
		case 23:		NV_Ith_S(abstol__, 23) = value_new;    break;
		case 24:		NV_Ith_S(abstol__, 24) = value_new;    break;
		case 25:		NV_Ith_S(abstol__, 25) = value_new;    break;
		case 26:		NV_Ith_S(abstol__, 26) = value_new;    break;
		case 27:		NV_Ith_S(abstol__, 27) = value_new;    break;
		case 28:		NV_Ith_S(abstol__, 28) = value_new;    break;
		case 29:		NV_Ith_S(abstol__, 29) = value_new;    break;
		case 30:		NV_Ith_S(abstol__, 30) = value_new;    break;
		case 31:		NV_Ith_S(abstol__, 31) = value_new;    break;
		case 32:		NV_Ith_S(abstol__, 32) = value_new;    break;
		case 33:		NV_Ith_S(abstol__, 33) = value_new;    break;
		case 34:		NV_Ith_S(abstol__, 34) = value_new;    break;
		case 35:		NV_Ith_S(abstol__, 35) = value_new;    break;
		case 36:		NV_Ith_S(abstol__, 36) = value_new;    break;
		case 37:		NV_Ith_S(abstol__, 37) = value_new;    break;
		case 38:		NV_Ith_S(abstol__, 38) = value_new;    break;
		case 39:		NV_Ith_S(abstol__, 39) = value_new;    break;
		case 40:		NV_Ith_S(abstol__, 40) = value_new;    break;
		default: fprintf(stderr,"ERROR - setAbstol - index = %d out of bounds\n",index);    break;
		}
	}

	double Solveode::calc_Bi(){
		return (pow((1.+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cai_old_),2.))),(-1.)));	//5
	}
	double Solveode::calc_Bss(){
		return (pow((1.+((CMDN_tot*Km_CMDN)/pow((Km_CMDN+Cass_old_),2.))),(-1.)));	//6
	}
	double Solveode::calc_BJSR(){
		return (pow((1.+((CSQN_tot*Km_CSQN)/pow((Km_CSQN+CaJSR_old_),2.))),(-1.)));	//7
	}
	double Solveode::calc_J_rel(){
		return ((v1*(P_O1_old_+P_O2_old_)*(CaJSR_old_-Cass_old_)*P_RyR_old_));	//8
	}
	double Solveode::calc_J_tr(){
		return (((CaNSR_old_-CaJSR_old_)/tau_tr));	//9
	}
	double Solveode::calc_J_xfer(){
		return (((Cass_old_-Cai_old_)/tau_xfer));	//10
	}
	double Solveode::calc_J_leak(){
		return ((v2*(CaNSR_old_-Cai_old_)));	//11
	}
	double Solveode::calc_J_up(){
		return (((v3*pow(Cai_old_,2.))/(pow(Km_up,2.)+pow(Cai_old_,2.))));	//12
	}
	double Solveode::calc_J_trpn(){
		return ((((k_plus_htrpn*Cai_old_*(HTRPN_tot-HTRPN_Ca_old_))+(k_plus_ltrpn*Cai_old_*(LTRPN_tot-LTRPN_Ca_old_)))-((k_minus_htrpn*HTRPN_Ca_old_)+(k_minus_ltrpn*LTRPN_Ca_old_))));	//13
	}
	double Solveode::calc_P_C1(){
		return ((1.-(P_C2_old_+P_O1_old_+P_O2_old_)));	//18
	}
	double Solveode::calc_i_CaL(){
		return ((g_CaL*O_old_*(V_old_-E_CaL)));	//21
	}
	double Solveode::calc_C1(){
		//return ((1-(O_old_+C2_old_+C2_old_+C3_old_+C4_old_+I1_old_+I2_old_+I3_old_)));ERRO	//23
                return ((1.0-(O_old_+C2_old_+C3_old_+C4_old_+I1_old_+I2_old_+I3_old_))); //23
	}
	double Solveode::calc_alpha(){
		return (((0.4*exp(((V_old_+12.)/10.))*((1.+(0.7*exp(((-pow((V_old_+40.),2.))/10.))))-(0.75*exp(((-pow((V_old_+20.),2.))/400.)))))/(1.+(0.12*exp(((V_old_+12.)/10.))))));	//30
	}
	double Solveode::calc_beta(){
		return ((0.05*exp(((-(V_old_+12.))/13.))));	//31
	}
	double Solveode::calc_gamma(){
		return (((Kpc_max*Cass_old_)/(Kpc_half+Cass_old_)));	//32
	}
	double Solveode::calc_Kpcf(){
		return ((13.*(1.-exp(((-pow((V_old_+14.5),2.))/100.)))));	//33
	}
	double Solveode::calc_i_pCa(){
		return (((i_pCa_max*pow(Cai_old_,2.))/(pow(Km_pCa,2.)+pow(Cai_old_,2.))));	//34
	}
	double Solveode::calc_i_NaCa(){
		//return ((((((((k_NaCa*1)/(pow(K_mNa,3)+pow(Nao,3)))*1)/(K_mCa+Cao))*1)/(1+(k_sat*exp((((eta-1)*V_old_*F)/(R*T))))))*((exp(((eta*V_old_*F)/(R*T)))*pow(Nai_old_,3)*Cao)-(exp((((eta-1)*V_old_*F)/(R*T)))*pow(Nao,3)*Cai_old_))));	//35
               double u1,u2,u3,u4; //implementacao dorothy
               u1=1./(pow(K_mNa,3.)+pow(Nao,3.));
               u2=1./(K_mCa+Cao);
               u3=1./(1.+k_sat*exp(((eta-1.)*V_old_*F)/(R*T)));
               u4=(exp((eta*V_old_*F)/(R*T))*pow(Nai_old_,3.)*Cao)-(exp((eta-1.)*V_old_*F/(R*T))*pow(Nao,3.)*Cai_old_);
               return (k_NaCa*u1*u2*u3*u4);

	}
	double Solveode::calc_i_Cab(){
		return ((g_Cab*(V_old_-calc_E_CaN())));	//36
	}
	double Solveode::calc_E_CaN(){
		return ((((R*T)/(2.*F))*log((Cao/Cai_old_))));	//37
	}
	double Solveode::calc_i_Na(){
		return ((g_Na*O_Na_old_*(V_old_-calc_E_Na())));	//39
	}
	double Solveode::calc_E_Na(){
		return ((((R*T)/F)*log((((0.9*Nao)+(0.1*Ko))/((0.9*Nai_old_)+(0.1*Ki_old_))))));	//40
	}
	double Solveode::calc_C_Na3(){
		return ((1.-(O_Na_old_+C_Na1_old_+C_Na2_old_+IF_Na_old_+I1_Na_old_+I2_Na_old_+IC_Na2_old_+IC_Na3_old_)));	//41
	}
	double Solveode::calc_alpha_Na11(){
		return ((3.802/((0.1027*exp(((-(V_old_+2.5))/17.)))+(0.2*exp(((-(V_old_+2.5))/150.))))));	//50
	}
	double Solveode::calc_alpha_Na12(){
		return ((3.802/((0.1027*exp(((-(V_old_+2.5))/15.)))+(0.23*exp(((-(V_old_+2.5))/150.))))));	//51
	}
	double Solveode::calc_alpha_Na13(){
		return ((3.802/((0.1027*exp(((-(V_old_+2.5))/12.)))+(0.25*exp(((-(V_old_+2.5))/150.))))));	//52
	}
	double Solveode::calc_beta_Na11(){
		return ((0.1917*exp(((-(V_old_+2.5))/20.3))));	//53
	}
	double Solveode::calc_beta_Na12(){
		return ((0.2*exp(((-(V_old_-2.5))/20.3))));	//54
	}
	double Solveode::calc_beta_Na13(){
		return ((0.22*exp(((-(V_old_-7.5))/20.3))));	//55
	}
	double Solveode::calc_alpha_Na3(){
		return ((7e-7*exp(((-(V_old_+7.))/7.7))));	//56
	}
	double Solveode::calc_beta_Na3(){
		return ((0.0084+(0.00002*(V_old_+7.))));	//57
	}
	double Solveode::calc_alpha_Na2(){
		return ((1./((0.188495*exp(((-(V_old_+7.))/16.6)))+0.393956)));	//58
	}
	double Solveode::calc_beta_Na2(){
		return (((calc_alpha_Na13()*calc_alpha_Na2()*calc_alpha_Na3())/(calc_beta_Na13()*calc_beta_Na3())));	//59
	}
	double Solveode::calc_alpha_Na4(){
		return ((calc_alpha_Na2()/1000.));	//60
	}
	double Solveode::calc_beta_Na4(){
		return (calc_alpha_Na3());	//61
	}
	double Solveode::calc_alpha_Na5(){
		return ((calc_alpha_Na2()/95000.));	//62
	}
	double Solveode::calc_beta_Na5(){
		return ((calc_alpha_Na3()/50.));	//63
	}
	double Solveode::calc_i_Nab(){
		return ((g_Nab*(V_old_-calc_E_Na())));	//64
	}
	double Solveode::calc_i_Kto_f(){
		return ((g_Kto_f*pow(ato_f_old_,3.)*ito_f_old_*(V_old_-calc_E_K())));	//66
	}
	double Solveode::calc_E_K(){
		return ((((R*T)/F)*log((Ko/Ki_old_))));	//67
	}
	double Solveode::calc_alpha_a(){
		return ((0.18064*exp((0.03577*(V_old_+30.)))));	//70
	}
	double Solveode::calc_beta_a(){
		return ((0.3956*exp(((-0.06237)*(V_old_+30.)))));	//71
	}
	double Solveode::calc_alpha_i(){
		return (((0.000152*exp(((-(V_old_+13.5))/7.)))/((0.0067083*exp(((-(V_old_+33.5))/7.)))+1.)));	//72
	}
	double Solveode::calc_beta_i(){
		return (((0.00095*exp(((V_old_+33.5)/7.)))/((0.051335*exp(((V_old_+33.5)/7.)))+1.)));	//73
	}
	double Solveode::calc_i_Kto_s(){
		return ((g_Kto_s*ato_s_old_*ito_s_old_*(V_old_-calc_E_K())));	//74
	}
	double Solveode::calc_ass(){
		return ((1./(1.+exp(((-(V_old_+22.5))/7.7)))));	//77
	}
	double Solveode::calc_iss(){
		return ((1./(1.+exp(((V_old_+45.2)/5.7)))));	//78
	}
	double Solveode::calc_tau_ta_s(){
		return (((0.493*exp(((-0.0629)*V_old_)))+2.058));	//79
	}
	double Solveode::calc_tau_ti_s(){
		return ((270.+(1050./(1.+exp(((V_old_+45.2)/5.7))))));	//80
	}
	double Solveode::calc_i_K1(){
		return (((((0.2938*Ko)/(Ko+210.))*(V_old_-calc_E_K()))/(1.+exp((0.0896*(V_old_-calc_E_K()))))));	//81
	}
	double Solveode::calc_i_Ks(){
		return ((g_Ks*pow(nKs_old_,2.)*(V_old_-calc_E_K())));	//82
	}
	double Solveode::calc_alpha_n(){
		//return ((0.00000481333*(V_old_+26.5)*(1-exp(((-0.128)*(V_old_+26.5))))));ERRO	//84
                return ((0.00000481333*(V_old_+26.5))/(1.0-exp(((-0.128)*(V_old_+26.5)))));      //31
	}
	double Solveode::calc_beta_n(){
		return ((0.0000953333*exp(((-0.038)*(V_old_+26.5)))));	//85
	}
	double Solveode::calc_i_Kur(){
		return ((g_Kur*aur_old_*iur_old_*(V_old_-calc_E_K())));	//86
	}
	double Solveode::calc_tau_aur(){
		return (((0.493*exp(((-0.0629)*V_old_)))+2.058));	//89
	}
	double Solveode::calc_tau_iur(){
		return ((1200.-(170./(1.+exp(((V_old_+45.2)/5.7))))));	//90
	}
	double Solveode::calc_i_Kss(){
		return ((g_Kss*aKss_old_*iKss_old_*(V_old_-calc_E_K())));	//91
	}
	double Solveode::calc_tau_Kss(){
		return (((39.3*exp(((-0.0862)*V_old_)))+13.17));	//94
	}
	double Solveode::calc_i_Kr(){
		return ((g_Kr*O_K_old_*(V_old_-(((R*T)/F)*log((((0.98*Ko)+(0.02*Nao))/((0.98*Ki_old_)+(0.02*Nai_old_))))))));	//95
	}
	double Solveode::calc_C_K0(){
		return ((1.-(C_K1_old_+C_K2_old_+O_K_old_+I_K_old_)));	//96
	}
	double Solveode::calc_alpha_a0(){
		return ((0.022348*exp((0.01176*V_old_))));	//101
	}
	double Solveode::calc_beta_a0(){
		return ((0.047002*exp(((-0.0631)*V_old_))));	//102
	}
	double Solveode::calc_alpha_a1(){
		return ((0.013733*exp((0.038198*V_old_))));	//103
	}
	double Solveode::calc_beta_a1(){
		return ((0.0000689*exp(((-0.04178)*V_old_))));	//104
	}
	double Solveode::calc_alpha_i_duplicated_rapid_delayed_rectifier_potassium_current(){
		return ((0.090821*exp((0.023391*(V_old_+5.)))));	//105
	}
	double Solveode::calc_beta_i_duplicated_rapid_delayed_rectifier_potassium_current(){
		return ((0.006497*exp(((-0.03268)*(V_old_+5.)))));	//106
	}
	double Solveode::calc_i_NaK(){
		return (((((i_NaK_max*calc_f_NaK()*1.)/(1.+pow((Km_Nai/Nai_old_),1.5)))*Ko)/(Ko+Km_Ko)));	//107
	}
	double Solveode::calc_f_NaK(){
		return ((1./(1.+(0.1245*exp((((-0.1)*V_old_*F)/(R*T))))+(0.0365*calc_sigma()*exp((((-V_old_)*F)/(R*T)))))));	//108
	}
	double Solveode::calc_sigma(){
		return (((1./7.)*(exp((Nao/67300.))-1.)));	//109
	}
	double Solveode::calc_i_ClCa(){
		return ((((g_ClCa*calc_O_ClCa()*Cai_old_)/(Cai_old_+Km_Cl))*(V_old_-E_Cl)));	//110
	}
	double Solveode::calc_O_ClCa(){
		return ((0.2/(1.+exp(((-(V_old_-46.7))/7.8)))));	//111
	}
	double Solveode::calc_Istim(){
                switch(protocolo) {
                     case 0:  return (ifnumber_0()); break;
                     case 1:  return (ifnumber_1()); break;
                     case 2:  return (ifnumber_2()); break;
                }
	}
	double Solveode::ifnumber_0(){
             //pacing
             if(((time_new>=IstimStart)&&(time_new<=IstimEnd)&&(((time_new-IstimStart)-(floor(((time_new-IstimStart)/IstimPeriod))*IstimPeriod))<=IstimPulseDuration))){
                     return (IstimAmplitude);
             }else{
                     return (0.0);
             }

            //voltage clamp simples
            /*double t1, t2, v_hold, v_volt;
            t1 = 50.0;
            t2 = 3500.0;
            v_hold = (-80.   - V_old_)/0.001;
            v_volt = (v_step - V_old_)/0.001;
            return -(v_hold*((time_new < t1) | (time_new > t2)) + v_volt*((time_new >= t1) & (time_new <= t2)));
            */
        }

        double Solveode::ifnumber_1(){
            //  voltage clamp pré-pulso - IK_total
            double t1, t2, t3, v_hold1, v_hold2, v_volt;
            t1 = 50.0;
            t2 = 100.0; 
            t3 = 3100.0;
            v_hold1 = (-80.0 - V_old_)/0.001;
            v_hold2 = (-40.0 - V_old_)/0.001;
            v_volt  = (v_step - V_old_)/0.001;
            return -(v_hold1*((time_new < t1) | (time_new > t3)) + v_hold2*((time_new >= t1) & (time_new <= t2)) + v_volt*((time_new > t2) & (time_new <= t3)));

/*            double t1, t2, v_hold, v_volt;
            t1 = 10.0;
            t2 = 5000.0;
            v_hold = (-80.0 - V_old_)/0.001;
            v_volt  = (v_step - V_old_)/0.001;
            if(time_new < t1)
               return v_hold;
            else if((time_new >= t1) && (time_new <= t2))
               return v_volt;*/
            //return -(v_hold*((time_new < t1)) + v_volt*((time_new >= t1) & (time_new <= t2)));
        }

        double Solveode::ifnumber_2(){
            //  Simple voltage clamp pré-pulso - ICa_total
            double t1, t2, t3, v_hold1, v_hold2, v_volt;
            t1 = 50.0;
            t2 = 100.0; 
            t3 = 350.0;
            v_hold1 = (-80.0 - V_old_)/0.001;
            v_hold2 = (-40.0 - V_old_)/0.001;
            v_volt  = (v_step - V_old_)/0.001;
            return -(v_hold1*((time_new < t1) | (time_new > t3)) + v_hold2*((time_new >= t1) & (time_new <= t2)) + v_volt*((time_new > t2) & (time_new <= t3)));

//             double t1, t2, v_hold, v_volt;
//             t1 = 10.0;
//             t2 = 250.0;
//             v_hold = (-80.0 - V_old_)/0.001;
//             v_volt  = (v_step - V_old_)/0.001;
// 
//             if(time_new < t1)
//                return v_hold;
//             else if((time_new >= t1) && (time_new <= t2))
//                return v_volt;
// 

            //return -(v_hold*((time_new < t1)) + v_volt*((time_new >= t1) & (time_new <= t2)));

        }


static int check_flag(void *flagvalue, char *funcname, int opt){
	int *errflag;
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1);}
	else if (opt == 1) {
		errflag = (int *) flagvalue;
		if (*errflag < 0) {
			fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",funcname, *errflag);
			return(1); }}
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
		return(1); }
	return 0;
}

static int f__(realtype time, N_Vector dependent_variable__, N_Vector dep_var_dot__, void *f_data__){
	Solveode *ode = (Solveode *) f_data__;
	ode->setVariables( 0 ,NV_Ith_S(dependent_variable__, 0));
	ode->setVariables( 1 ,NV_Ith_S(dependent_variable__, 1));
	ode->setVariables( 2 ,NV_Ith_S(dependent_variable__, 2));
	ode->setVariables( 3 ,NV_Ith_S(dependent_variable__, 3));
	ode->setVariables( 4 ,NV_Ith_S(dependent_variable__, 4));
	ode->setVariables( 5 ,NV_Ith_S(dependent_variable__, 5));
	ode->setVariables( 6 ,NV_Ith_S(dependent_variable__, 6));
	ode->setVariables( 7 ,NV_Ith_S(dependent_variable__, 7));
	ode->setVariables( 8 ,NV_Ith_S(dependent_variable__, 8));
	ode->setVariables( 9 ,NV_Ith_S(dependent_variable__, 9));
	ode->setVariables( 10 ,NV_Ith_S(dependent_variable__, 10));
	ode->setVariables( 11 ,NV_Ith_S(dependent_variable__, 11));
	ode->setVariables( 12 ,NV_Ith_S(dependent_variable__, 12));
	ode->setVariables( 13 ,NV_Ith_S(dependent_variable__, 13));
	ode->setVariables( 14 ,NV_Ith_S(dependent_variable__, 14));
	ode->setVariables( 15 ,NV_Ith_S(dependent_variable__, 15));
	ode->setVariables( 16 ,NV_Ith_S(dependent_variable__, 16));
	ode->setVariables( 17 ,NV_Ith_S(dependent_variable__, 17));
	ode->setVariables( 18 ,NV_Ith_S(dependent_variable__, 18));
	ode->setVariables( 19 ,NV_Ith_S(dependent_variable__, 19));
	ode->setVariables( 20 ,NV_Ith_S(dependent_variable__, 20));
	ode->setVariables( 21 ,NV_Ith_S(dependent_variable__, 21));
	ode->setVariables( 22 ,NV_Ith_S(dependent_variable__, 22));
	ode->setVariables( 23 ,NV_Ith_S(dependent_variable__, 23));
	ode->setVariables( 24 ,NV_Ith_S(dependent_variable__, 24));
	ode->setVariables( 25 ,NV_Ith_S(dependent_variable__, 25));
	ode->setVariables( 26 ,NV_Ith_S(dependent_variable__, 26));
	ode->setVariables( 27 ,NV_Ith_S(dependent_variable__, 27));
	ode->setVariables( 28 ,NV_Ith_S(dependent_variable__, 28));
	ode->setVariables( 29 ,NV_Ith_S(dependent_variable__, 29));
	ode->setVariables( 30 ,NV_Ith_S(dependent_variable__, 30));
	ode->setVariables( 31 ,NV_Ith_S(dependent_variable__, 31));
	ode->setVariables( 32 ,NV_Ith_S(dependent_variable__, 32));
	ode->setVariables( 33 ,NV_Ith_S(dependent_variable__, 33));
	ode->setVariables( 34 ,NV_Ith_S(dependent_variable__, 34));
	ode->setVariables( 35 ,NV_Ith_S(dependent_variable__, 35));
	ode->setVariables( 36 ,NV_Ith_S(dependent_variable__, 36));
	ode->setVariables( 37 ,NV_Ith_S(dependent_variable__, 37));
	ode->setVariables( 38 ,NV_Ith_S(dependent_variable__, 38));
	ode->setVariables( 39 ,NV_Ith_S(dependent_variable__, 39));
	ode->setVariables( 40 ,NV_Ith_S(dependent_variable__, 40));
	ode->setParameters(0,time);
	double *t = ode->solveDiff();
	NV_Ith_S(dep_var_dot__, 0) = t[0];
	NV_Ith_S(dep_var_dot__, 1) = t[1];
	NV_Ith_S(dep_var_dot__, 2) = t[2];
	NV_Ith_S(dep_var_dot__, 3) = t[3];
	NV_Ith_S(dep_var_dot__, 4) = t[4];
	NV_Ith_S(dep_var_dot__, 5) = t[5];
	NV_Ith_S(dep_var_dot__, 6) = t[6];
	NV_Ith_S(dep_var_dot__, 7) = t[7];
	NV_Ith_S(dep_var_dot__, 8) = t[8];
	NV_Ith_S(dep_var_dot__, 9) = t[9];
	NV_Ith_S(dep_var_dot__, 10) = t[10];
	NV_Ith_S(dep_var_dot__, 11) = t[11];
	NV_Ith_S(dep_var_dot__, 12) = t[12];
	NV_Ith_S(dep_var_dot__, 13) = t[13];
	NV_Ith_S(dep_var_dot__, 14) = t[14];
	NV_Ith_S(dep_var_dot__, 15) = t[15];
	NV_Ith_S(dep_var_dot__, 16) = t[16];
	NV_Ith_S(dep_var_dot__, 17) = t[17];
	NV_Ith_S(dep_var_dot__, 18) = t[18];
	NV_Ith_S(dep_var_dot__, 19) = t[19];
	NV_Ith_S(dep_var_dot__, 20) = t[20];
	NV_Ith_S(dep_var_dot__, 21) = t[21];
	NV_Ith_S(dep_var_dot__, 22) = t[22];
	NV_Ith_S(dep_var_dot__, 23) = t[23];
	NV_Ith_S(dep_var_dot__, 24) = t[24];
	NV_Ith_S(dep_var_dot__, 25) = t[25];
	NV_Ith_S(dep_var_dot__, 26) = t[26];
	NV_Ith_S(dep_var_dot__, 27) = t[27];
	NV_Ith_S(dep_var_dot__, 28) = t[28];
	NV_Ith_S(dep_var_dot__, 29) = t[29];
	NV_Ith_S(dep_var_dot__, 30) = t[30];
	NV_Ith_S(dep_var_dot__, 31) = t[31];
	NV_Ith_S(dep_var_dot__, 32) = t[32];
	NV_Ith_S(dep_var_dot__, 33) = t[33];
	NV_Ith_S(dep_var_dot__, 34) = t[34];
	NV_Ith_S(dep_var_dot__, 35) = t[35];
	NV_Ith_S(dep_var_dot__, 36) = t[36];
	NV_Ith_S(dep_var_dot__, 37) = t[37];
	NV_Ith_S(dep_var_dot__, 38) = t[38];
	NV_Ith_S(dep_var_dot__, 39) = t[39];
	NV_Ith_S(dep_var_dot__, 40) = t[40];
	return 0;
}



float __agos_factorial(int f){
	if(f>=0 & f<2)
		return 1.0;
	else if(f < 0)
		return 0.0/0.0;
	for(int i=f-1; i>=2; i--)
		f *= i;
	return (float)f;
}
