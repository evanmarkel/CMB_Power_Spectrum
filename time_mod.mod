	  ¸*  o   k820309    Ö          15.0        kÖqU                                                                                                           
       time_mod.f90 TIME_MOD                                                    
                                                          
                                                          
                                                          
                                                                                                                                                                                                                                                                    
                 
                    JxÞ±A        2.99792458D8                                                 
                   
                  )ÿmrìD<                                                    	     
                 
                 ðk$	Â?        8.3D-5#         @                                  
                	   #ODEINT%SIZE    #ODEINT%SIGN    #ODEINT%ABS    #YSTART    #X1    #X2    #EPS    #H1    #HMIN    #DERIVS    #RKQS    #OUTPUT %                                                    SIZE                                                  SIGN                                                  ABS           
                                                  
               &                                                     
                                      
                
                                      
                
                                      
                
                                      
                
                                      
      #         @                                        	               #X    #Y    #DYDX              
                                     
                
                                                   
              &                                                                                                       
               &                                           #         @                                        	            	   #Y    #DYDX    #X    #HTRY    #EPS    #YSCAL    #HDID    #HNEXT     #DERIVS !             
                                                  
               &                                                     
                                                   
 	             &                                                     
                                    
                 
                                     
                
                                     
                
                                                   
 
             &                                                                                         
                                                      
       #         @                                   !     	               #X "   #Y #   #DYDX $             
                                "     
                
                                #                   
              &                                                                                    $                   
               &                                           #         @                                   %     	               #X &   #Y '             
                                &     
                
                                '                   
              &                                           #         @    @                             (                	   #BSSTEP%SIZE )   #BSSTEP%MAXVAL *   #BSSTEP%ABS +   #BSSTEP%MAX ,   #BSSTEP%MIN -   #Y .   #DYDX /   #X 0   #HTRY 1   #EPS 2   #YSCAL 3   #HDID 4   #HNEXT 5   #DERIVS 6                                               )     SIZE                                             *     MAXVAL                                             +     ABS                                             ,     MAX                                             -     MIN           
                               .                   
 
              &                                                     
                                 /                   
              &                                                     
                                0     
                 
                                 1     
                
                                 2     
                
                                 3                   
              &                                                                                     4     
                                                 5     
       #         @                                   6     	               #X 7   #Y 8   #DYDX 9             
                                7     
                
                                8                   
              &                                                                                    9                   
 	              &                                           #         @                                  :                   #SPLINE%SIZE ;   #X <   #Y =   #YP1 >   #YPN ?   #Y2 @                                               ;     SIZE           
                                <                   
              &                                                     
                                 =                   
              &                                                     
                                 >     
                
                                 ?     
                                                @                   
               &                                                                                       A     
                 
                 VM¡ #D        3.08568025D22                                            B     
                 
                 ¹?        .1D0                                            C     
                   
                  ¬<kM$Ë:                                                    D     
                 
                 Ù?        .4D0                                            E     
                 
                                 0.D0%         @                               F                   
       #SPLINT%SIZE G   #SPLINT%MAX H   #SPLINT%MIN I   #XA J   #YA K   #Y2A L   #X M                                               G     SIZE                                             H     MAX                                             I     MIN           
                                J                   
              &                                                     
                                 K                   
              &                                                     
                                 L                   
              &                                                     
                                 M     
                 @                                N                     @ @                              O                   
                &                                                    @                                P                   
                &                                                    @                                Q                   
                &                                                      @                                R                     @ @                              S                   
                &                                                    @ @                              T                   
                &                                                    @ @                              U                   
                &                                           #         @                                   V                    #INITIALIZE_TIME_MOD%SQRT W   #INITIALIZE_TIME_MOD%EXP X   #INITIALIZE_TIME_MOD%LOG Y                                             W     SQRT                                           X     EXP                                           Y     LOG #         @    @                             Z                    #X [   #Y \   #DETADA ]             
  @                              [     
                
                                 \                   
 	             &                                                     D                                ]                   
 
              &                                           #         @    @                             ^                    #X _   #Y `             
                                 _     
                
                                 `                   
              &                                           %         @                               a                   
       #GET_H%SQRT b   #GET_H%EXP c   #X d                                             b     SQRT                                           c     EXP           
  @                              d     
      %         @                               e                   
       #GET_H_P%SQRT f   #GET_H_P%EXP g   #X h                                             f     SQRT                                           g     EXP           
  @                              h     
      %         @                               i                   
       #GET_DH_P%SQRT j   #GET_DH_P%EXP k   #X l                                             j     SQRT                                           k     EXP           
  @                              l     
      %         @                               m                    
       #X_IN n             
  @                              n     
                   fn#fn    ¾   @   J   HEALPIX_TYPES    þ   @   J   PARAMS    >  @   J   SPLINE_1D_MOD    ~  @   J   ODE_SOLVER "   ¾  p       I4B+HEALPIX_TYPES !   .  p       DP+HEALPIX_TYPES      |       C+PARAMS      p       H_0+PARAMS      v       OMEGA_R+PARAMS "      Ó       ODEINT+ODE_SOLVER '   Ó  =      ODEINT%SIZE+ODE_SOLVER '     =      ODEINT%SIGN+ODE_SOLVER &   M  <      ODEINT%ABS+ODE_SOLVER )        a   ODEINT%YSTART+ODE_SOLVER %     @   a   ODEINT%X1+ODE_SOLVER %   U  @   a   ODEINT%X2+ODE_SOLVER &     @   a   ODEINT%EPS+ODE_SOLVER %   Õ  @   a   ODEINT%H1+ODE_SOLVER '     @   a   ODEINT%HMIN+ODE_SOLVER )   U  `      ODEINT%DERIVS+ODE_SOLVER +   µ  @   a   ODEINT%DERIVS%X+ODE_SOLVER +   õ     a   ODEINT%DERIVS%Y+ODE_SOLVER .        a   ODEINT%DERIVS%DYDX+ODE_SOLVER '   	        ODEINT%RKQS+ODE_SOLVER )   ¬	     a   ODEINT%RKQS%Y+ODE_SOLVER ,   8
     a   ODEINT%RKQS%DYDX+ODE_SOLVER )   Ä
  @   a   ODEINT%RKQS%X+ODE_SOLVER ,     @   a   ODEINT%RKQS%HTRY+ODE_SOLVER +   D  @   a   ODEINT%RKQS%EPS+ODE_SOLVER -        a   ODEINT%RKQS%YSCAL+ODE_SOLVER ,     @   a   ODEINT%RKQS%HDID+ODE_SOLVER -   P  @   a   ODEINT%RKQS%HNEXT+ODE_SOLVER .     `      ODEINT%RKQS%DERIVS+ODE_SOLVER 0   ð  @   a   ODEINT%RKQS%DERIVS%X+ODE_SOLVER 0   0     a   ODEINT%RKQS%DERIVS%Y+ODE_SOLVER 3   ¼     a   ODEINT%RKQS%DERIVS%DYDX+ODE_SOLVER )   H  V      ODEINT%OUTPUT+ODE_SOLVER +     @   a   ODEINT%OUTPUT%X+ODE_SOLVER +   Þ     a   ODEINT%OUTPUT%Y+ODE_SOLVER    j  ó       BSSTEP+BS_MOD #   ]  =      BSSTEP%SIZE+BS_MOD %     ?      BSSTEP%MAXVAL+BS_MOD "   Ù  <      BSSTEP%ABS+BS_MOD "     <      BSSTEP%MAX+BS_MOD "   Q  <      BSSTEP%MIN+BS_MOD          a   BSSTEP%Y+BS_MOD #        a   BSSTEP%DYDX+BS_MOD     ¥  @   a   BSSTEP%X+BS_MOD #   å  @   a   BSSTEP%HTRY+BS_MOD "   %  @   a   BSSTEP%EPS+BS_MOD $   e     a   BSSTEP%YSCAL+BS_MOD #   ñ  @   a   BSSTEP%HDID+BS_MOD $   1  @   a   BSSTEP%HNEXT+BS_MOD %   q  `      BSSTEP%DERIVS+BS_MOD '   Ñ  @   a   BSSTEP%DERIVS%X+BS_MOD '        a   BSSTEP%DERIVS%Y+BS_MOD *        a   BSSTEP%DERIVS%DYDX+BS_MOD %   )         SPLINE+SPLINE_1D_MOD *   ª  =      SPLINE%SIZE+SPLINE_1D_MOD '   ç     a   SPLINE%X+SPLINE_1D_MOD '   s     a   SPLINE%Y+SPLINE_1D_MOD )   ÿ  @   a   SPLINE%YP1+SPLINE_1D_MOD )   ?  @   a   SPLINE%YPN+SPLINE_1D_MOD (        a   SPLINE%Y2+SPLINE_1D_MOD      }       MPC+PARAMS      t       OMEGA_B+PARAMS    ü  p       RHO_C+PARAMS    l  t       OMEGA_M+PARAMS $   à  t       OMEGA_LAMBDA+PARAMS %   T  ¡       SPLINT+SPLINE_1D_MOD *   õ  =      SPLINT%SIZE+SPLINE_1D_MOD )   2  <      SPLINT%MAX+SPLINE_1D_MOD )   n  <      SPLINT%MIN+SPLINE_1D_MOD (   ª     a   SPLINT%XA+SPLINE_1D_MOD (   6     a   SPLINT%YA+SPLINE_1D_MOD )   Â     a   SPLINT%Y2A+SPLINE_1D_MOD '   N  @   a   SPLINT%X+SPLINE_1D_MOD      @       N_T    Î         X_T    Z         A_T    æ         Z_T    r   @       N_ETA    ²          X_ETA    >!         ETA    Ê!         ETA2 $   V"          INITIALIZE_TIME_MOD )   ö"  =      INITIALIZE_TIME_MOD%SQRT (   3#  <      INITIALIZE_TIME_MOD%EXP (   o#  <      INITIALIZE_TIME_MOD%LOG    «#  b       DERIVS_ETA    $  @   a   DERIVS_ETA%X    M$     a   DERIVS_ETA%Y "   Ù$     a   DERIVS_ETA%DETADA    e%  V       OUTPUT_ETA    »%  @   a   OUTPUT_ETA%X    û%     a   OUTPUT_ETA%Y    &  v       GET_H    ý&  =      GET_H%SQRT    :'  <      GET_H%EXP    v'  @   a   GET_H%X    ¶'  z       GET_H_P    0(  =      GET_H_P%SQRT    m(  <      GET_H_P%EXP    ©(  @   a   GET_H_P%X    é(  |       GET_DH_P    e)  =      GET_DH_P%SQRT    ¢)  <      GET_DH_P%EXP    Þ)  @   a   GET_DH_P%X    *  Z       GET_ETA    x*  @   a   GET_ETA%X_IN 