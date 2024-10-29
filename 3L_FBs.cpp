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
		static bool   qg = 0, ql = 0, qs = 0, qll1 = 0, qlr1 = 0, qll2 = 0, qlr2 = 0; // Switchs
		static bool   q[7] = {0,0,0,0,0,0,0}; // Switch Array (Dimension = NºSwitch)
		
		// Inputs (S denotes save variables, read)
		static double egs = 0., igs = 0., vCs = 0., vtris = 0.;
	
		// Scalar PWM
		static double vm1 = 0., vM1 = 0., vm2 = 0., vM2 = 0., vm3 = 0., vM3 = 0., vM_min = 0., vm_max = 0., u = 0.;
		static double vu = 0., vumax = 0., vumin = 0.;			// Apportionment factor
		static double vg0_ref = 0., vl0_ref = 0., vs0_ref = 0.; // 3L
		static double vll10_ref = 0., vlr10_ref = 0., vll20_ref = 0., vlr20_ref = 0.; // FBs
		
		// Bus Control - Simple PI 
		static float  vCs_ref = 170.; 							// Reference Bus Voltage
		static double vCs_error = 0.; 							// Error Bus Voltage
		static double Ig_ref = 0.;						    	// Reference Current Amplitude 
		static double Ivc_error = 0., Pvc_error = 0., Ig0 = 3.; // PI Error
		static float  kpv = 0.2, kiv = 20.; 					// Controller Gains
		
		// Grid Current Control - Resonant PI
		static double ig_ref = 0., vg_ref = 0.;					// Reference
		static double ig_error = 0., ig_error_p = 0.;			// PIR Error

		static double thetasw_ref = 0.;							//  Period-switched phase
		static double F11s = 0., F12s = 0., F21s = 0., F22s = 0., H11s = 0., H21s = 0.;
		static double X1ccs = 0., X1cs = 0., X2ccs = 0., X2cs = 0.; 
		static double kis = 1000., kps = 10;				    // Controller Gains	
		
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
			vCs_error = vCs_ref - vCs;					// Bus Voltage Error
			
			Ivc_error = Ivc_error + kiv*hpwm*vCs_error; // Integral Error
			Pvc_error = kpv*vCs_error;                  // Proportional Error
			
			Ig_ref = Pvc_error + Ivc_error + Ig0; 		// PI Out - Grid Current Amplitude
			
			if(Ig_ref > 30) Ig_ref = 30;
			if(Ig_ref < -30) Ig_ref = -30;
			
			// ######### Bus Control - End #########
			
			
			
			// ######### Grid Current Control - Start ######### 
			ig_ref = Ig_ref*sin(theta);
			ig_error = (igs - ig_ref);
			
			// Cut Frequency
			kis = 1000.;//150.;1000.;
			kps = 10.;//1.;10.;
			
			// Constants
			wg = 2.0*pi*fg;
			thetasw_ref = wg*hpwm;
			F11s = cos(thetasw_ref); 
			F12s = sin(thetasw_ref)/wg;
			F21s = -wg*sin(thetasw_ref);
			F22s = F11s;
			H11s =  2.0*kis*sin(thetasw_ref)/wg;
			H21s = (cos(thetasw_ref)-1.0)*2.0*kis;
			
			// Controller 
			X1ccs = X1cs;  
			X2ccs = X2cs;  
			X1cs = F11s*X1ccs + F12s*X2ccs + H11s*ig_error;
			X2cs = F21s*X1ccs + F22s*X2ccs + H21s*ig_error;
			vg_ref = X1cs + ig_error*kps ;  // Output PI
			//vg_ref = 109.5461*sqrt(2)*sin(theta-0.1598);
			
			// Saturator
			if (vg_ref >= 1.5*vCs_ref)
		    {
		        X1cs = 1.5*vCs_ref;
		        vg_ref = 1.5*vCs_ref;
		    }
		    else if (vg_ref <= -1.5*vCs_ref)
		    {
	            X1cs = -1.5*vCs_ref;
	            vg_ref = -1.5*vCs_ref;
		    }
		    
			// ######### Grid Current Control - End ######### 
			
			
			////////////////////
			// il signal assignment through powers to avoid phase shift above 60°
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
			vl_ref = vlL + Rl*ils + Ll*(ils_ref - ils)*fs;
//			vl_ref = 110.*sqrt(2)*sin(theta); // Open Loop
			
			
			// Saturator
			if (vl_ref >= 1.5*vCs_ref) vl_ref = 1.5*vCs_ref;
		    else if (vl_ref <= -1.5*vCs_ref) vl_ref = -1.5*vCs_ref;
			
			// ######### Load Current Control - End #########			
			
			
			
			///////////////////// PWM ////////////////////////
			// 3L PWM
			// Scalar Limits
			vM1 =  vg_ref;
			vM2 =  vl_ref;
			vM3 =  0;
			
			vm1 =  vg_ref;
			vm2 =  vl_ref;
			vm3 =  0;
			
			//min(vM1,vM2,vM3)
			if((vM1<=vM2)&&(vM1<=vM3)) vM_min = vM1;
			else{if(vM2<=vM3) vM_min = vM2; else vM_min = vM3;}
			
			//max(vm1,vm2,vm3)
			if((vm1>=vm2)&&(vm1>=vm3)) vm_max = vm1;
			else{if(vm2>=vm3) vm_max = vm2; else vm_max = vm3;}	 
			
			// Scalar PWM
			vumax =  vCs_ref/2. - vm_max; 
			vumin = -vCs_ref/2. - vM_min;
			u = 0.5;
			vu = u*vumin + (1-u)*vumax; 
			//vu = 0;
			vs0_ref = (vu)/vCs_ref + 0.5;
			vl0_ref = (vu + vl_ref)/vCs_ref + 0.5;
			vg0_ref = (vu + vg_ref)/vCs_ref + 0.5;
			
		}
		
		if(vg0_ref   >= vtris) q[0] = 1.0; else q[0] = 0.0;
		if(vl0_ref   >= vtris) q[1] = 1.0; else q[1] = 0.0;
		if(vs0_ref   >= vtris) q[2] = 1.0; else q[2] = 0.0;	
		if(vll10_ref >= vtris) q[3] = 1.0; else q[3] = 0.0;
		if(vlr10_ref >= vtris) q[4] = 1.0; else q[4] = 0.0;
		if(vll20_ref >= vtris) q[5] = 1.0; else q[5] = 0.0;
		if(vlr20_ref >= vtris) q[6] = 1.0; else q[6] = 0.0;
				
		qg   = q[0];
	    ql   = q[1];
		qs   = q[2];
		qll1 = q[3];
		qlr1 = q[4];
		qll2 = q[5];
		qlr2 = q[6];
		
//		qll1 = 1.0;
//		qlr1 = 0.0;
//		qll2 = 1.0;
//		qlr2 = 0.0;

		
		///////////////////////////////////
		// Output
		
		out[0]  =  qg;
		out[1]  =  ql;
		out[2]  =  qs;
		out[3]  =  qll1;
		out[4]  =  qlr1;
		out[5]  =  qll2;
		out[6]  =  qlr2;
		
		out[7]  =  vg0_ref;
		out[8]  =  vl0_ref;
		out[9]  =  vs0_ref;
		//out[7]  =  Il_ref;
		//out[8]  =  sign;
		//out[9]  =  Ig_ref;
		
		out[10] =  vll10_ref;
		out[11] =  vlr10_ref;
		out[12] =  vll20_ref;
		out[13] =  vlr20_ref;
		out[14] =  vg_ref;
		out[15] =  vl_ref;
		out[16] =  ig_ref;
		out[17] =  v2L1_ref;
		out[18] =  v2L2_ref;
		out[19] =  V2L1_ref;
		
		*pnError = 0; //Success
		// ............end	
	}
}

