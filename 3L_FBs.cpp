#include "dll.h"
#include <windows.h>
#include <math.h>

#define hpwm 1E-04
#define N 167


extern "C"
{
		 
	/////////////////////////////////////////////////////////////////////
	// FUNCTION: SimulationStep
	//   This function runs at every time step.
	//double t: (read only) time
	//double delt: (read only) time step as in Simulation control
	//double *in: (read only) zero based array of input values. in[0] is the first node, in[1] second input...
	//double *out: (write only) zero based array of output values. out[0] is the first node, out[1] second output...
	//int *pnError: (write only)  assign  *pnError = 1;  if there is an error and set the error message in szErrorMsg
	//    strcpy(szErrorMsg, "Error message here..."); 
	DLLIMPORT void SimulationStep(
			double t, double delt, double *in, double *out,
			 int *pnError, char * szErrorMsg,
			 void ** ptrUserData, int nThreadIndex, void * pAppPtr)
	{
		// Constants
		static float  fg = 60., wg = 0., pi = 3.1415, fs = 1/hpwm;
		
		// Processing Loop Control
		static float  td = 0.;
	
		// Switch Status
		static bool   qse = 0, qsh = 0, qll1 = 0, qlr1 = 0, qll2 = 0, qlr2 = 0; // Switchs
		static bool   q[6] = {0,0,0,0,0,0}; // Switch Array (Dimension = N?Switch)
		
		// Inputs (S denotes save variables, read)
		static double egs = 0., igs = 0., vCs = 0., vtris = 0.;
	
		// Scalar PWM
		static double vse_ref = 0., vsh_ref = 0.,vse0_ref = 0., vsh0_ref = 0.; // 3L
		static double vll10_ref = 0., vlr10_ref = 0., vll20_ref = 0., vlr20_ref = 0.; // FBs

		// Controle do Barramento
		float kpigs = 0.5, kiigs = 20.0;                // Ganhos proporcional e integral
		float e_vdcs, vCs_ref=510;                      // Erro de tensão, referência e tensão medida no barramento
		float Vcts, Vctfs = 0.0;                        // Valor temporário e valor de controle com função integral
		float Ig_refs, ig_refs;                         // Tempo de amostragem, corrente de referência e corrente ajustado

		// Controle da Corrente da Rede
		static double erroigs = 0.0, erro1cs = 0.0;               // Erros de corrente (atual e anterior)
		static double varefs = 0.0, vg_refs = 0.0;                // Saída do controlador e referência de tensão

		static double kis = 1000.0, kps = 10.0;                   // Ganhos do controlador PI ressonante
		static double ws = 0.0, f = 0.0, tetas = 0.0;             // Frequência angular, frequência e ajuste de fase
		static double F11s = 0.0, F12s = 0.0, F21s = 0.0, F22s = 0.0; // Coeficientes do controlador
		static double H11s = 0.0, H21s = 0.0;                     // Coeficientes de ganho integral do controlador

		static double X1ccs = 0.0, X2ccs = 0.0;                   // Estados anteriores
		static double X1cs = 0.0, X2cs = 0.0;                     // Estados atuais
		
		// Load  Current Control - Predictive
		static double vl_ref = 0., vlL = 0.;
		static double Rl = 0.1, Ll = 5e-3; // Inductance Parameters
		
		// Phase Detection Control - PLL
		static double theta = 0., wf = 0, Pd = 0.0;
		static double X1t = 0.0, dX1t = 0.0;
		static double X1tf = 0.0, Pdf = 0.0, e_pd = 0.0;
		static double X1a = 0., X2a = 0.0, dX1a = 0.0, dX2a = 0.0, X1af = 0.0, X2af = 0.0, kix = 0.0, kpx = 0.0;
		static double A1 = -8.796*5., A2 = -39.48*5.*5., A3 = 1.0, A4 = 0.0, B1 = 1., B2= 0.0, C1 = 0.0, C2 = 39.48*5.*5., D1 = 0.0;
 	
		// Grid Voltage RMS Detector
		static int EG[N] = {}, N_pos = 0;
		static double ceg = 0., ceg_rms = 0., ceg_amplitude = 0.;

		// FBs Bus Control - FBs Voltage - Simple PI
		static double vC2L1s = 0., vC2L2s = 0.;
		static double vC2L1s_ref = 100., vC2L2s_ref = 70.;
		static double vC2L_ref = (vC2L1s_ref+vC2L2s_ref)/2;
		static double vC2L1s_error = 0., vC2L2s_error = 0.;
		static double v2L1_ref = 0., v2L2_ref = 0., V2L1_ref = 0., V2L2_ref = 0.;
		static double kpvl = 0.1, kivl = 20.; // Same gains for both PIs
		static double Iv2L1_error = 0, Iv2L2_error = 0, Pv2L1_error = 0, Pv2L2_error = 0; // Errors integral and proportional
		static double P2L1 = 0, P2L2 = 0, Pltot = 0; // 2L Bus Power
		static double iC2L1s = 0, iC2L2s = 0; // 2L Bus Current
		
		// FBs Bus Control - Load Current - Simple PI
		static double ils = 0., ils_ref = 0.;
		static double Il_ref = 15, Il_ref0 = 5; // IL* = Eg*/S* = 9.1 A = 1 p.u; S* = 1000 W  --- 13.6 -> 1.5 p.u.
		static double v2Ls_max = 0., v2Ls_ref_max = 0., v2Ls_error_max = 0.;
		static double kpil = 1.0, kiil = 10; // Same gains for both PIs
		static double Iil_error = 0, Pil_error = 0; // Errors integral and proportional
		static int sign = 0;
		
		
		
		///////////////////////////////////
		// Input
		egs 	     = in[0];	// Grid Voltage
		igs          = in[1];	// Grid Current
		vCs		  	 = in[2];   // 3L Bus Voltage
		vtris  		 = in[3];	// Triangular
		ils 		 = in[4];	// load Current
		vC2L1s 	 	 = in[5];   // 2L1 Bus Voltage
		vC2L2s 		 = in[6];   // 2L2 Bus Voltage
		iC2L1s       = in[7];	// 2L1 Bus Current
		iC2L2s       = in[8];	// 2L2 Bus Current
		
		//if(deg==1) Il_ref = 26.9;
		
		if(td <= t) // Start Process
		{
			td = td + hpwm;
			
		 	///////////////////////PLL - Grid Voltage - Phase Detection////////////////////////////
			kix = 1.15;
			kpx = 0.15;
		
			Pd = egs*cos(theta);
		    
			// Low pass filter
			X1a = X1af;
			X2a = X2af;
			
			dX1a = A1*X1a + A2*X2a + B1*Pd;
			dX2a = A3*X1a + A4*X2a + B2*Pd;
			
			X1af = X1a + dX1a*hpwm;
			X2af = X2a + dX2a*hpwm;
			
			Pdf = C1* X1af + C2*X2af + D1*Pd; // Filter output
			
			e_pd = -(0. - Pdf); // Error
		
			// PI Controller + Integrator			  
			X1t = X1tf;
			X1tf = X1t + kix*e_pd*hpwm;
			dX1t = X1tf + kpx*e_pd;
			wf = 2*pi*fg + dX1t;
			
			theta = theta + wf*hpwm; // PLL output
			   
			if(theta>=2.*pi) theta -= 2.*pi; 		// Restart Phase
			if(theta<= 0.) theta += 2.*pi;			// Restart Phase
			
			
			
			////////////////////// Cumulative Grid RMS Voltage ////////////////////////////
			ceg -= EG[N_pos]*EG[N_pos];
			ceg += egs*egs;
			
			EG[N_pos] = egs;         		    // Instantly saves grid voltage
			N_pos++;							// Voltage array progress index
			
			ceg_rms = sqrt(ceg/N);				// Grid RMS Voltage
			ceg_amplitude = ceg_rms*sqrt(2);    // Grid Peak Voltage
			
			if(N_pos >= N) N_pos -= N;				// Restart Position
			
			
			
			///////////////////// Control Sequence ////////////////////////
			// ######### Bus Control - Start ######### 
			e_vdcs = vCs_ref - vCs;
			
			Vcts = Vctfs;
			
			Vctfs = Vcts + kiigs*e_vdcs*hpwm;
			
			if(Vctfs > 40.) Vctfs = 40.;
			if(Vctfs < -40.) Vctfs = -40.;			
			
			Ig_refs = Vctfs + kpigs*e_vdcs;
			
			if(Ig_refs > 40.) Ig_refs = 40.;
			if(Ig_refs < 0.) Ig_refs = 0.;
			
			ig_refs = Ig_refs*sin(theta);
			// ######### Bus Control - End #########


			// ######### Grid Current Control - Start ######### 
			erroigs = (-ig_refs + igs);
			
			//Frequencia de Corte do controlador
			kis = 1000.;//150 100 1000;
			kps = 10.;//10 1;
			
			//Constantes do Cotrolador
			f = fg;
			ws = 2.0*pi*f;
			tetas = ws*hpwm;
			F11s = cos(tetas); 
			F12s = sin(tetas)/ws;
			F21s = -ws*sin(tetas);
			F22s = F11s;
			H11s =  2.0*kis*sin(tetas)/ws;
			H21s = (cos(tetas)-1.0)*2.0*kis;
			
			//CONTROLADORES DAS CORRENTES ig1 
			
			X1ccs = X1cs;  
			X2ccs = X2cs;  
			X1cs = F11s*X1ccs + F12s*X2ccs + H11s*erro1cs;
			X2cs = F21s*X1ccs + F22s*X2ccs + H21s*erro1cs;
			varefs = X1cs + erroigs*kps;  // Saida 

			/**/
			if (varefs >= 1.5*vCs_ref)
		    {
		        X1cs = 1.5*vCs_ref;
		        varefs = 1.5*vCs_ref;
		    }
		    else
		    {
		        if (varefs <= -1.5*vCs_ref)
		        {
		            X1cs = -1.5*vCs_ref;
		            varefs = -1.5*vCs_ref;
		        }
		    }
		    
		    
			erro1cs = erroigs; //Erro anterior
			vsh_ref = varefs;
			// ######### Grid Current Control - End ######### 
			
			////////////////////
			// il signal assignment through powers to avoid phase shift above 60?
			P2L1 = vC2L1s_ref*iC2L1s;
			P2L2 = vC2L2s_ref*iC2L2s;
			Pltot = P2L1+P2L2;
			
			if(Pltot > 0) sign = 1;
			else sign = -1;
			////////////////////
			
			// ######### FBs Bus Voltage Control + FB PWM - Begin #########
			// Simple PI
			//////////////////// 2L1
			vC2L1s_error = sign*(vC2L1s_ref - vC2L1s);
//			vC2L1s_error = (vC2L1s_ref - vC2L1s);
			Iv2L1_error = Iv2L1_error + kivl*hpwm*vC2L1s_error;
			Pv2L1_error = kpvl*vC2L1s_error;
			V2L1_ref = (Iv2L1_error + Pv2L1_error);
			
			v2L1_ref = V2L1_ref*sin(theta);
			
			//////////////////// 2L2
			vC2L2s_error = sign*(vC2L2s_ref - vC2L2s);
//			vC2L2s_error = (vC2L2s_ref - vC2L2s);
			Iv2L2_error = Iv2L2_error + kivl*hpwm*vC2L2s_error;
			Pv2L2_error = kpvl*vC2L2s_error;
			V2L2_ref = (Iv2L2_error + Pv2L2_error);
			
			v2L2_ref = V2L2_ref*sin(theta);
			
			//////////////////// FB PWM
			// FB1 PWM
			vll10_ref =  (v2L1_ref/2)/vC2L1s_ref + 0.5;
			vlr10_ref = -(v2L1_ref/2)/vC2L1s_ref + 0.5;
			
			// FB2 PWM	
			vll20_ref =  (v2L2_ref/2)/vC2L2s_ref + 0.5;
			vlr20_ref = -(v2L2_ref/2)/vC2L2s_ref + 0.5;			

			// ######### FBs Bus Voltage Control + FB PWM - End #########
			
			
			// ######### Load Current Control - Start ######### 			
			if(abs(V2L1_ref) > abs(V2L2_ref))
			{
				v2Ls_max = V2L1_ref;
				v2Ls_ref_max = 0.9*vC2L1s_ref; // Modulation Index  ma = 0.9
			}
			else 
			{
				v2Ls_max = V2L2_ref;
				v2Ls_ref_max = 0.9*vC2L2s_ref; // Modulation Index  ma = 0.9
			}
			
			// PI
			v2Ls_error_max = -(abs(v2Ls_ref_max) - abs(v2Ls_max));
			Iil_error = Iil_error + kiil*hpwm*v2Ls_error_max;
			Pil_error = kpil*v2Ls_error_max;
			Il_ref = sign*(Iil_error + Pil_error + Il_ref0);
//			Il_ref = (Iil_error + Pil_error + Il_ref0);
			
			// By the Power
//			if(abs(iC2L1s) > abs(iC2L2s)) //P2L1/vC2L1s_ref = iC2L1s
//				Il_ref = sqrt(2)*abs(iC2L1s);
//			else				//iC2L1s < iC2L2s
//				Il_ref = sqrt(2)*abs(iC2L2s);
				
//			Saturator
			if(sign > 0)
			{
				if(Il_ref >= 50) Il_ref = 50;
    			else if(Il_ref <= 0) Il_ref = 0;
			}
			else
			{
				if(Il_ref <= -50) Il_ref = -50;
		    	else if(Il_ref >= 0) Il_ref = 0;
			}
			
			// Predictive Control 
			ils_ref = Il_ref*sin(theta); // Resistive Load (fp = 1)
			vlL = v2L1_ref + v2L2_ref;
			//vl_ref = vlL + Rl*ils + Ll*(ils_ref - ils)*fs;
			vl_ref = 110.*sqrt(2)*sin(theta); // Open Loop
			vse_ref = vl_ref - egs;
			
			
			// Saturator
			if (vl_ref >= 1.5*vCs_ref) vl_ref = 1.5*vCs_ref;
		    else if (vl_ref <= -1.5*vCs_ref) vl_ref = -1.5*vCs_ref;
			
			// ######### Load Current Control - End #########			
			
			
			
			///////////////////// PWM ////////////////////////
			// 3L PWM
			vsh0_ref = (-vsh_ref)/vCs_ref + 0.5;
			vse0_ref = (-vse_ref)/vCs_ref + 0.5;
			
			
		}
		
		
		if(vse0_ref   >= vtris) q[0] = 1.0; else q[0] = 0.0;
		if(vsh0_ref   >= vtris) q[1] = 1.0; else q[1] = 0.0;	
		if(vll10_ref >= vtris) q[2] = 1.0; else q[2] = 0.0;
		if(vlr10_ref >= vtris) q[3] = 1.0; else q[3] = 0.0;
		if(vll20_ref >= vtris) q[4] = 1.0; else q[4] = 0.0;
		if(vlr20_ref >= vtris) q[5] = 1.0; else q[5] = 0.0;
				
		
	    qse   = q[0];
		qsh   = q[1];
		qll1 = q[2];
		qlr1 = q[3];
		qll2 = q[4];
		qlr2 = q[5];
		
//		qll1 = 1.0;
//		qlr1 = 0.0;
//		qll2 = 1.0;
//		qlr2 = 0.0;

		
		///////////////////////////////////
		// Output
		
		
		out[1]  =  qse;
		out[2]  =  qsh;
		out[3]  =  qll1;
		out[4]  =  qlr1;
		out[5]  =  qll2;
		out[6]  =  qlr2;
		
		//out[7]  =  vg0_ref;
		//out[8]  =  vl0_ref;
		//out[9]  =  vs0_ref;
		//out[7]  =  Il_ref;
		//out[8]  =  sign;
		//out[9]  =  Ig_ref;
		
		out[10] =  vll10_ref;
		out[11] =  vlr10_ref;
		out[12] =  vll20_ref;
		out[13] =  vlr20_ref;
		out[14] =  vsh_ref;
		out[15] =  vse_ref;
		out[16] =  ig_refs;
		out[17] =  v2L1_ref;
		out[18] =  v2L2_ref;
		out[19] =  V2L1_ref;
		
		*pnError = 0; //Success
		// ............end	
	}
}

