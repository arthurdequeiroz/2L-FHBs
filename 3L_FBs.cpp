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
		//............begin
		
		static double f = 60., wff = 0., fc = 10000.;
		static double pi = 3.1415;
		static double td = 0.;
		static double vsh_ref = 0., vse_ref = 0., Vc1 = 0., ig = 0., il = 0.;
	
	
		static double qg = 0., ql = 0., qs = 0., qx = 1; 
		static int q[5] = {0,0,0,0,0};
		static double vC = 0.;
	
		static double vm1=0.0,vM1=0.0,vm2=0.0,vM2=0.0,vm3=0.0,vM3=0.0,vM_min=0.0,vm_max=0., migt = 0., vx = 0., vs0 = 0, vl0 = 0.;
		
	
	
		static double ig_ref = 0., theta = 0.;
		static double wf = 0.; 
		
		static double A1 = -8.796*5., A2 = -39.48*5.*5., A3 = 1.0, A4 = 0.0;
		static double B1 = 1., B2= 0.0, C1 = 0.0, C2 = 39.48*5.*5., D1 = 0.0;
		static double Pd = 0.0;
		
		static double X1t = 0.0, dX1t = 0.0;
		static double X1tf = 0.0, Pdf = 0.0, e_pd = 0.0;
		
		static double X1a = 0., X2a = 0.0, dX1a = 0.0, dX2a = 0.0;
		static double X1af = 0.0, X2af = 0.0, kix = 0.0, kpx = 0.0;
		
		static double ws = 0.;
		
		static double erroigs = 0., kis = 0., kps = 0., tetas = 0.;
		
		static double F11s = 0., F12s = 0., F21s = 0., F22s = 0., H11s = 0., H21s = 0.;
		static double varefs = 0.;
	
		static double e_vdcs = 0., Vcts = 0., Vctfs = 0., kiigs = 0., kpigs = 0.;
		static double vCs = 0, Ig_refs = 0., ig_refs=0., ig1s =0, X1ccs = 0., X1cs = 0.,X2ccs = 0., X2cs = 0., erro1cs =0, vCpos=0,vCneg=0; 
		static double vg_refs = 0., vg0 = 0.;
	
		static double vxmax = 0., vxmin = 0., Ts = 0., vbase = 0, trian = 0., Rs = 0., Ls = 0.;
		
		static double x = 0, flag = 1., flagx = 0, flag_swell = 0, trigger = 0;
		static double vC_ref = 510, wffL= 0;
	
		static double iss = 0, is_ref = 0, is_error = 0, Il, il_ref;
	
		static double egs = 0, igs = 0, Vgm = 0, Vgp = 0, Vgm_ref = 0, Tm = 0;
		static double El_ref = 220, VG[N] = {};
		
		static int va = 0;
		static float vg_soma = 0;
		static float vg_rms = 0;
		
		static float t0 = 0;
		
//		static float td2 = 0, derv = 0;
		
		Ts = hpwm;
		wff = 2.*pi*60.;
		wffL = 2.*pi*60.;
		vbase = 110*sqrt(2);
	
		///////////////////////////////////
		// ENTRADAS
		
		trian    = in[0];
		egs	     = in[1];
		igs  	 = in[2];
		vCs      = in[3]; 
		

		
		if(td <= t)
		{
			
			td = td + hpwm;
			
				
		 	///////////////////////PLL EG////////////////////////////
			
			kix = 1.15;
			kpx = 0.15;
		
			Pd = (egs/2)*cos(theta);
		    //FPB
			X1a = X1af;
			X2a = X2af;
			
			dX1a = A1*X1a + A2*X2a + B1*Pd;
			dX2a = A3*X1a + A4*X2a + B2*Pd;
			
			X1af = X1a + dX1a*hpwm;
			X2af = X2a + dX2a*hpwm;
			
			Pdf = C1* X1af + C2*X2af + D1*Pd; //Saida do Filtro
			//Erro
			e_pd = -(0. - Pdf);
		
			//Controlador PI + Integrador			  
			X1t = X1tf;
			X1tf = X1t + kix*e_pd*hpwm;
			dX1t = X1tf + kpx*e_pd;
			wf = wff + dX1t;
			
			theta = theta + wf*hpwm; //Saida da PLL
			   
			if(theta>=2.*pi) theta = theta - 2.*pi;
			if(theta<= 0.) theta = theta + 2.*pi;
						
		
			//kpigs = 0.2;//0.05;//0.2;
			kpigs = 0.5;//0.2
			kiigs = 20;//20;
		    
			e_vdcs = vC_ref - vCs;
			
			Vcts = Vctfs;
			
			Vctfs = Vcts + kiigs*e_vdcs*hpwm;
			
			if(Vctfs > 40.) Vctfs = 40.;
			if(Vctfs < -40.) Vctfs = -40.;			
			
			Ig_refs = Vctfs + kpigs*e_vdcs;
			
			if(Ig_refs > 40.) Ig_refs = 40.;
			if(Ig_refs < 0.) Ig_refs = 0.;
			
			ig_refs = Ig_refs*sin(theta);
			//Reguladores de corrente (PI sin) - Ig1 
			
			//erro1s = (-ig1_refs + is);
			
			// PR RESSONANTE COMEÇA AQUI
			erroigs = (-ig_refs + igs);
			
			//Frequencia de Corte do controlador
			kis = 1000.;//150 100 1000;
			kps = 10.;//10 1;
			
			//Constantes do Cotrolador
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
			if (varefs >= 1.5*vC_ref)
		    {
		        X1cs = 1.5*vC_ref;
		        varefs = 1.5*vC_ref;
		    }
		    else
		    {
		        if (varefs <= -1.5*vC_ref)
		        {
		            X1cs = -1.5*vC_ref;
		            varefs = -1.5*vC_ref;
		        }
		    }
		    
		    
			erro1cs = erroigs; //Erro anterior
			vg_refs = varefs;
			
			/*   TENSÕES DE REFERÊNCIA   */		


			vsh_ref = vg_refs;

			vse_ref = (egs - 110.*sqrt(2)*sin(theta));
			
			
			vl0 = (-vse_ref);			
			vs0 = (-vsh_ref);			

		}
		
		/*   PWM - DEFINIÇÃO DO ESTADO DAS CHAVES   */		
		if(vl0 >= trian) q[0] = 1.0; else q[0] = 0.0;
		if(vs0 >= trian) q[1] = 1.0; else q[1] = 0.0;	
		
		ql  = q[0];
		qs  = q[1];
		
		/*   SAÍDAS   */
		out[0] =  qs;
		out[1] =  ql;
		out[2] =  ig_refs; // Flag
		out[3] =  vsh_ref;
	
		
		*pnError = 0; //Success
		// ............end	
	}
}














