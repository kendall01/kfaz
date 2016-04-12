function [ Ma,T,T0 ] = findMaTempsVar( mdot,A,Tm,P0,P,RF )
% ASSUME: Variable Cp,Cv
% FIND: Ma, T & T0
% NOTE: Approach from ME140 Lecture 3 page 21,22

error = 0.01;

% (a) Assume: T0=Tm initially
T0guess = Tm;
i = 1;
j = 1;

while i==1;
    % (b) Find MFP via Eq. (1)
    MFP = (mdot/A)*sqrt(R*T0guess)/P0;

    % (c) Determine T, P0/P, and Ma that satisfies Eqn. (2),(3),&(4)    
        % (c1) Assume Tguess = Tm initially
        Tguess = Tm;
        
        while j ==1;
            [Cp,Cv,k] = sp_heats(Tguess);
            % (c2) Evaluate 3 integrals: I1,I2,& cp,ave
            % ---- from T0guess to Tguess
            I1 = integral_1(T0guess,Tguess);
            I2 = integral_2(T0guess,Tguess);
            cp_ave = integral_3(T0guess,Tguess);
        
            % (c3) Find P0/P via Eq. (2)
            P0byP = exp(I2);
        
            % (c4) Find Ma that satisfies Eq. (3)
            Ma = solve( MFP == ( M*sqrt(k)/P0byP )*( 1+ ((k*R)/(2*cp_ave))*M^2 )^0.5, M );
        
            % (c5) Find T via Eq. (4), if T = Tguess then T=answer, 
            % ---- else let Tguess = T & repeat from (c2).
            T = Tm/( 1+ RF*Ma^2*( (k*R)/(2*cp_ave) ) );
            
            if abs(T-Tguess) < error
                j = 0;     
            else
                Tguess = T;
            end
        end
        
    % (d) Find T0 via Eq. (5), if T0=T0guess T0=answer , else T0guess = T0
    % --- and repeat from (b).
    T0 = T*( 1+ Ma^2*( (k*R)/(2*cp_ave) ) );

    if abs(T0-T0guess)<error
        i = 0; 
    else
        T0guess = T0;
    end
    

end
% --- Else, let T0guess = T0 and repeat from (b)
        T0guess = T0;
end

