!     Last change:  JCH  23 Mar 2012    3:26 pm_v1
!ENDOWMENT ECONOMY. COST OF DEFAULTING: 1 PERIOD EXCLUSION AND UTILITY LOSS. SHOCK TO RISK PREMIA
!IF indicator market access = 1 ==> GOOD, LOW RISK PREMIUM
!IF indicator market access = 2 ==> BAD, HIGH RISK PREMIUM
module param
!INCLUDE 'link_fnl_shared.h'
DOUBLE PRECISION :: beta, sigma, output_cost, r, b_inf, b_sup, b_inf_short, b_sup_short, y_inf, y_sup, eps_inf, eps_sup, std_eps,&
                    rho, mean_y, std_y, pi_number, zero, width, prob_excl_end, pi, coupon, delta, gamma, d0, d1, cdf_inf, &
                    cdf_sup, phi, escala, r_long, omega, psi, tax_inf, tax_sup, net_inf, net_sup, b_recov,&
                    kappa, disut_default, alpha, absortion
                    
         

parameter (std_eps = 0.0340d+0, rho = 0.66d+0, mean_y = -0.5d+0*std_eps**2,  pi_number = 3.1415926535897932d+0, &
           zero = 0d+0,  width = 1.5d+0, cdf_inf = 3.167d-5, cdf_sup = 1d+0 - 3.167d-5, pi = 0.00d+0,  &
          gamma = 0d+0, r = 0.04, phi = 0.0d+0, escala = 1d+2, kappa = 0.0d+0, &
          alpha = 0.667d+0, absortion = 0.12d+0)

INTEGER :: b_num_short, b_num_long, y_num, eps_num, quad_num, tax_num, access_num, nout, i_y_global, i_excl_global, cdf_num, &
           i_b_short_global, i_b_long_global, i_default_global, i_b_short_next, i_b_long_next,&
           quad_num_v, quad_num_q, b_num_short_finer, b_num_long_finer, y_num_finer, i_access_global, i_access_initial,&
           num_short_next, y_num_coarse, b_num_long_coarse, b_num_short_coarse, net_num
parameter (b_num_short = 40, b_num_long = 40, y_num =30, quad_num_v = 100, access_num = 2, eps_num = quad_num_v, quad_num_q = 100, cdf_num =50,&
           b_num_short_finer = 100, b_num_long_finer = 100, y_num_finer =30, num_short_next = 50, tax_num = 10, &
           y_num_coarse =15, b_num_long_coarse = 30, b_num_short_coarse = 30, net_num = 100)

DOUBLE PRECISION :: b_grid_long(1:b_num_long), b_grid_short(1:b_num_short), g_grid(access_num), y_grid(1:y_num), default_grid(1:2), indicator_tirar,&
                    vector_v_global(y_num), y_initial, g_initial, b_short_initial, b_long_initial, b_global,counter, b_long_global, &
                    b_short_global, q_short_global, cdf_grid(cdf_num), quad_w_v(1:quad_num_v), quad_x_v(1:quad_num_v),&
                    quad_w_hermite(1:eps_num), quad_x_hermite(1:eps_num), quad_w_q(1:quad_num_q), quad_x_q(1:quad_num_q),&
                    b_short_global_optimize, b_grid_long_finer(1:b_num_long_finer), b_grid_short_finer(1:b_num_short_finer), &
                    y_grid_finer(1:y_num_finer), b_next_short_grid(num_short_next), &
                    b_grid_long_coarse(1:b_num_long), b_grid_short_coarse(1:b_num_short), y_grid_coarse(1:y_num), indicator_global_search, cons_global, &
                    premia_grid(2), tax_opt_global, b_next_short_global, b_next_long_global, tax_grid(tax_num), c_global, g_global, output_global,&
                    net_grid(net_num), net_global, q_global, b_short_next_global, util_global, ev_global

DOUBLE PRECISION, DIMENSION(b_num_short, b_num_long, y_num, access_num) :: v_matrix
INTEGER, DIMENSION(b_num_short, b_num_long, y_num, 2) :: default_decision
DOUBLE PRECISION, DIMENSION(b_num_short, b_num_long, y_num, access_num,2) :: b_next_matrix
DOUBLE PRECISION, DIMENSION(b_num_short, b_num_long, y_num, access_num) :: v0_matrix, b0_next_short_matrix, b1_next_short_matrix,&
 q_paid_matrix, b0_next_long_matrix, v1_matrix, q_nodef_matrix, g0_matrix, g1_matrix, tax0_matrix, tax1_matrix
DOUBLE PRECISION, DIMENSION(2, 2) :: trans_matrix

DOUBLE PRECISION, DIMENSION(net_num, y_num_finer) :: tax0_interp_matrix, tax1_interp_matrix



!needed to use BS2IN - spline interpolation in 2 dimensions
INTEGER :: KORDER
PARAMETER( KORDER=3 )

INTEGER :: dim_auxiliar_long, dim_auxiliar_short, dim_auxiliar_long_finer, dim_auxiliar_short_finer
PARAMETER(dim_auxiliar_long = KORDER + b_num_long, dim_auxiliar_short = KORDER + b_num_short,&
          dim_auxiliar_long_finer = KORDER + b_num_long_finer, dim_auxiliar_short_finer = KORDER + b_num_short_finer )

DOUBLE PRECISION, DIMENSION(b_num_short_finer, y_num_finer, access_num) :: break_matrix_b_limit
DOUBLE PRECISION, DIMENSION(4, b_num_short_finer, y_num_finer, access_num) :: coeff_matrix_b_limit


DOUBLE PRECISION, DIMENSION (dim_auxiliar_long) :: b_long_knots
DOUBLE PRECISION, DIMENSION (dim_auxiliar_short) :: b_short_knots

DOUBLE PRECISION, DIMENSION (dim_auxiliar_long_finer) :: b_long_knots_finer
DOUBLE PRECISION, DIMENSION (dim_auxiliar_short_finer) :: b_short_knots_finer

DOUBLE PRECISION, DIMENSION(b_num_short, b_num_long, y_num, access_num) :: coeff_matrix_v0, coeff_matrix_q, coeff_matrix_v1
DOUBLE PRECISION, DIMENSION(b_num_short_finer, b_num_long_finer, y_num_finer, access_num) :: coeff_matrix_ev, coeff_matrix_q_menu, coeff_matrix_ev_excl


!DOUBLE PRECISION, DIMENSION(b_num_short, y_num, 2) :: break_matrix_v1
!DOUBLE PRECISION, DIMENSION(4, b_num_short, y_num, 2) :: coeff_matrix_v1

!DOUBLE PRECISION, DIMENSION(b_num_short_finer, y_num_finer, 2) :: break_matrix_ev_excl
!DOUBLE PRECISION, DIMENSION(4, b_num_short_finer, y_num_finer, 2) :: coeff_matrix_ev_excl

 end module




!SPECIFY GRID VALUES
subroutine compute_grid
USE param
DOUBLE PRECISION :: dos, DNORDF, delta_y, prob_mass, b_inf0, b_sup0, b_inf1, b_sup1, y_left, &
prob_vector(y_num), y_right, min_output, max_debt_min_revenue, tax_rate
INTEGER :: b_num_half, i,j
EXTERNAL :: DNORDF

coupon = (r + delta) / (1d+0 + r)
std_y = std_eps/ SQRT(1 - rho**2)

y_inf = mean_y - 5*std_y
y_sup = mean_y + 5*std_y
delta_y = 5*std_y / (y_num - 1d+0)

b_inf = 0.00001 !MIN(0.55d+0*EXP(y_sup), .99*EXP(y_inf))
b_sup =  min(1d+0, exp(y_inf)/coupon-0.05d+0)

b_inf_short = 0d+0 !MIN(0.55d+0*EXP(y_sup), .99*EXP(y_inf))
b_sup_short =  0.40d+0


open (10, FILE='b_grid_short.txt',STATUS='replace')
open (11, FILE='y_grid.txt',STATUS='replace')
open (12, FILE='b_grid_long.txt',STATUS='replace')
open (13, FILE='bounds.txt',STATUS='replace')


WRITE(13, '(F15.11, X, F15.11)') b_inf, b_sup
WRITE(13, '(F15.11, X, F15.11)') b_inf_short, b_sup_short
WRITE(13, '(F15.11, X, F15.11)') y_inf, y_sup


do i=1,b_num_short
   b_grid_short(i) = b_inf_short + (b_sup_short - b_inf_short)*(i-1d+0)/(b_num_short - 1d+0)
end do

do i=1,b_num_long
   b_grid_long(i) = b_inf + (b_sup - b_inf)*(i-1d+0)/(b_num_long - 1d+0)
end do


do i=1,b_num_short
 WRITE(10, '(F12.8, X, F12.8)') b_grid_short(i)
end do

do i=1,b_num_long
 WRITE(12, '(F12.8, X, F12.8)') b_grid_long(i)
end do

do i=1,y_num
   y_grid(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num - 1)
   WRITE(11, '(F12.8)') y_grid(i)
end do



close(10)
CLOSE(11)
CLOSE(12)
CLOSE(13)


open (16, FILE='trans_matrix.txt',STATUS='replace')
do i=1,access_num
   do j=1,access_num
       WRITE(16, '(F12.8)') trans_matrix(i,j)
   end do
end do
CLOSE(16)


do i=1,cdf_num
   cdf_grid(i) = cdf_inf + (cdf_sup - cdf_inf) * (i-1) / (cdf_num - 1)
end do

dos = 2d+00


default_grid(1) = 0.0d+0
default_grid(2) = 1.0d+0

eps_inf = -4*std_eps
eps_sup = 4*std_eps

open (20, FILE='b_short_finer.txt',STATUS='replace')
open (22, FILE='b_long_finer.txt',STATUS='replace')

do i=1,b_num_short_finer
   b_grid_short_finer(i) = b_inf_short + (b_sup_short - b_inf_short)*(i-1d+0)/(b_num_short_finer - 1d+0)
   WRITE(20, '(F12.8)') b_grid_short_finer(i)
end do

do i=1,b_num_long_finer
   b_grid_long_finer(i) = b_inf + (b_sup - b_inf)*(i-1d+0)/(b_num_long_finer - 1d+0)
   WRITE(22, '(F12.8)') b_grid_long_finer(i)
end do

close(20)
close(22)
do i=1,y_num_finer
   y_grid_finer(i) = y_inf + (y_sup - y_inf) * (i-1) / (y_num_finer - 1)
end do

!DEFIINE GRID OF RESERVE WHICH OPTIMIZE USES WHEN SOLVING THE OPTIMUM AT EACH RESERVE POINT
!IT IS DEFINED HERE SO THAT IT SIMPLIFIES THE FIRST LINES OF CODE IN OPTIMIZE.
!NEED TO CONSIDER POINTS CLOSE TO THE BOUNDARIES IN ORDER TO AVOID THE CODE TO PICK UP OPTIMAL
!VALUES AT THE BOUNDARIES WHEN THE OPTIMAL COMBINATION OF RESERVES, DEBT ARE AWAY FROM THE BOUNDARIES
!FOR RESERVES.
b_next_short_grid(1) = b_inf_short
b_next_short_grid(2) = b_inf_short + 0.001d+0
do i=2,num_short_next-3
   b_next_short_grid(i+1) = b_inf_short + (b_sup_short - b_inf_short) * (i-1)/ (num_short_next - 3)
end do
b_next_short_grid(num_short_next-1) = b_sup_short-0.001d+0
b_next_short_grid(num_short_next) = b_sup_short

open (16, FILE='pdf_y.txt',STATUS='replace')
!Aproximate unconditional p.d.f. of y
y_left  = (y_grid(1) + delta_y)/ std_y
y_right = (y_grid(y_num) - delta_y)/ std_y

prob_vector(1)     = DNORDF(y_left)
prob_vector(y_num) = 1d+0 - DNORDF(y_right)
write(16, '(F12.8)') prob_vector(1)
   do j=2,y_num-1
     y_right = (y_grid(j) + delta_y )/ std_y
     y_left  = (y_grid(j) - delta_y )/ std_y
     prob_vector(j) = (DNORDF(y_right) - DNORDF(y_left) ) !/prob_may_num
!     WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)')  y_grid(j), y_left, y_right, prob_vector(j)
     write(16, '(F12.8)') prob_vector(j)
   end do
write(16, '(F12.8)') prob_vector(y_num)
close (16)

end subroutine



!COMPUTE QUADRATURE POINTS AND WEIGHTS USING LEGENDRE AND GAUSSIAN QUADRATURE RULES
!THEY ARE USED TO COMPUTE NUMERICAL INTEGRALS
subroutine quadrature
USE param
INTEGER :: N_v, N_q, N, IWEIGH, IWEIGH4, NFIX
parameter(N_v = quad_num_v, N_q = quad_num_q, N = eps_num)
DOUBLE PRECISION:: ALFA, BETA1, QW_v(1:N_v), QX_v(1:N_v), QW_q(1:N_q), QX_q(1:N_q),QW(1:N), QX(1:N), QXFIX(2)
PARAMETER(ALFA=0, BETA1=0, IWEIGH=1, IWEIGH4=4)
external DGQRUL

!quad_w = weights using Gauss legendre quadrature rule
!quad_x = points using Gauss legendre quadrature rule
!quad_w_hermite = weights using Gauss Hermite quadrature rule
!quad_x_hermite = points using Gauss Hermite quadrature rule

NFIX = 0
CALL DGQRUL(N_v,IWEIGH, ALFA, BETA1, NFIX, QXFIX, QX_v, QW_v)

quad_w_v=QW_v
quad_x_v=QX_v

CALL DGQRUL(N_q,IWEIGH, ALFA, BETA1, NFIX, QXFIX, QX_q, QW_q)

quad_w_q=QW_q
quad_x_q=QX_q


CALL DGQRUL(N,IWEIGH4, ALFA, BETA1, NFIX, QXFIX, QX, QW)
quad_x_hermite = QX
quad_w_hermite = QW


end subroutine



DOUBLE PRECISION function dif_fun(y)
USE param
DOUBLE PRECISION, INTENT(IN) :: y
DOUBLE PRECISION :: v0_fun, v1_fun
EXTERNAL v0_fun, v1_fun


dif_fun = v0_fun(b_short_global, b_long_global,  y, i_access_global) - v1_fun(b_short_global, b_long_global, y, i_access_global)
!WRITE(nout, '(F12.8, X, F12.8, X, F12.8)') y, v0_fun(b_global, y), v1_fun(b_global, y)
end

!FUNCTION USED TO COMPUTE THE MINIMUM INCOME FOR WHICH THE CURRENT PARTY IN POWER DOES NOT DEFAULT
!b = outstanding debt
DOUBLE PRECISION function y_fun(b_short, b_long, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long
INTEGER, INTENT(IN) :: index_access
INTEGER :: MAXFN, num_tirar
PARAMETER(num_tirar = 100)
DOUBLE PRECISION :: dif_fun, ERRABS, ERRREL, left, right, y_max, y_min, dif_right, dif_left,&
                    y_tirar(num_tirar), tirar
EXTERNAL dif_fun, DZBREN

MAXFN=1000
ERRREL = 1D-10
ERRABS = 1D-10


y_max = rho*y_initial + (1-rho)*mean_y + width*eps_sup
y_min = rho*y_initial + (1-rho)*mean_y + width*eps_inf


b_short_global = b_short     !b_global IS USED IN dif_fun TO FIX THE VALUE OF b (OUTSTANDING DEBT)
b_long_global = b_long     !b_global IS USED IN dif_fun TO FIX THE VALUE OF b (OUTSTANDING DEBT)
i_access_global = index_access


dif_right = dif_fun(y_max)
dif_left  = dif_fun(y_min)

!1) DETERMINE WHETHER THERE IS AN INTERIOR ROOT OR NOT
!   a) IF THERE IS AN INTERIOR ROOT, SPECIFY WHETHER THE difference function IS INCREASING IN g OR NOT.
!   b) IF THERE IS NO INTERIOR ROOT, DETERMINE WHETHER THE GOV'T ALWAYS OR NEVER DEFAULTS

!WRITE(nout, *) y_min, y_max, dif_right, dif_left


if (dif_right*dif_left<0) then  !THERE IS AN INTERIOR ROOT
   left =  y_min
   right = y_max
   CALL DZBREN (dif_fun, ERRABS, ERRREL, left, right, MAXFN)
   y_fun = right
else if (dif_left>0) then
      y_fun = y_min
else !dif_right<0
      y_fun = y_max
end if

end



DOUBLE PRECISION function q_fun(b_short_next, b_long_next, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short_next, b_long_next, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION :: prob_no_def_tomorrow, y_current_type, DNORDF, y_fun, standarized_inf, standarized_sup,&
                    y_threshold, y_min, y_max, exp_q(2), scalar, y_next, y_next_def, q_paid_fun, &
                    DCSVAL, cdf, cdf_threshold, DNORIN, m_next, exp_m(2)

INTEGER :: i, NINTV, index_access_next
EXTERNAL DNORDF, y_fun, q_paid_fun, DCSVAL, DNORIN


exp_m = 0d+0
exp_q = 0d+0
scalar = 1d+0 / SQRT(pi_number)
NINTV = cdf_num - 1


if (b_long_next<=0) then
   q_fun = 100000 !coupon * (1d+0 - ((1d+0 - delta)/EXP(r))**counter)/ (delta + EXP(r)-1)
ELSEIF(b_long_next <=b_sup) THEN !SET POSITIVE PRICES ONLY FOR b' THAT ARE ABOVE THE MINIMUM VALUE.
                            !(AVOID OPTIMAL VALUES OUTSIDE THE RANGE)
   do index_access_next = 1,access_num

      y_threshold = y_fun(b_short_next, b_long_next, index_access_next)
      y_current_type = (y_threshold - rho * y - (1-rho)*mean_y) / std_eps
      cdf_threshold = DNORDF(y_current_type)

         if (cdf_threshold <= cdf_inf) then !THRESHOLD TOMORROW IS LOW
            prob_no_def_tomorrow = 1d+0
!            WRITE(nout, *) 'prob no def tomorrow = ', prob_no_def_tomorrow
            do i=1,quad_num_q
                cdf = 0.5*(quad_x_q(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
!                m_next = exp(-r) !1/(1+r)
                m_next = EXP(-r - premia_grid(index_access) * ( y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premia_grid(index_access)**2d+0)*std_eps**2d+0)
                exp_m(index_access_next) = exp_m(index_access_next) + m_next * coupon * 0.5* (cdf_sup - cdf_inf)*quad_w_q(i)
                exp_q(index_access_next) = exp_q(index_access_next) + (1d+0-delta) * m_next* &
                q_paid_fun(b_short_next, b_long_next, y_next, index_access_next) *0.5*(cdf_sup-cdf_inf)*quad_w_q(i)
                
            end do
            !NEED TO ADJUST EXPECTATION TO TAKE INTO ACCOUNT THAT THE
            !ACTUAL PROBABILITY MASS BETWEEN y_threshold and infty IS UNDERESTIMATED USING LEBESGE-..
      elseif (cdf_threshold >= cdf_sup) then !THRESHOLD TOMORROW IS HIGH
           exp_q(index_access_next) = 0d+0
      else
            prob_no_def_tomorrow = (1d+0 - cdf_threshold)
            do i=1,quad_num_q

                !compute integral for y > y_threshold
                cdf = 0.5*(quad_x_q(i) + 1)*(cdf_sup - cdf_threshold) + cdf_threshold
                y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                m_next = EXP(-r - premia_grid(index_access) * ( y_next - ( (1-rho)*mean_y + rho* y_initial)) - 0.5d+0*(premia_grid(index_access)**2d+0)*std_eps**2d+0)
                exp_m(index_access_next) = exp_m(index_access_next) +  m_next * coupon * 0.5* (cdf_sup - cdf_threshold)*quad_w_q(i)
        exp_q(index_access_next) = exp_q(index_access_next) + (1d+0 - delta) * m_next * &
        q_paid_fun(b_short_next, b_long_next, y_next, index_access_next) * 0.5* (cdf_sup - cdf_threshold)*quad_w_q(i)

            end do

      end if
      exp_q(index_access_next) = trans_matrix(index_access, index_access_next) * exp_q(index_access_next) / (cdf_sup - cdf_inf)
      exp_m(index_access_next) = trans_matrix(index_access, index_access_next) * exp_m(index_access_next) / (cdf_sup - cdf_inf)

   end do


q_fun = MAX(SUM(exp_m + exp_q), 0d+0)

end if


end


    


DOUBLE PRECISION function objective_function(X)
USE param
DOUBLE PRECISION, INTENT (IN) :: X(2)
INTEGER :: other_type, i,j,h, t, d, NINTV
DOUBLE PRECISION :: q_fun, acum, exp_v_next, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, q, DNORDF, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, cdf, cdf1,&
                    cdf_threshold, DNORIN, q_risk_free, output, b_short, b_long, q_short, EV_fun,&
                    q_menu_fun, penalty, indicator_constraint_violation, tax, tfp, labor, g, consum, &
                    util, interpolate, net_revenue, tirar, b_sup_local, b_limit_fun_spline, tax0_fun 

EXTERNAL q_fun, DNORDF, y_fun, v0_fun, v1_fun, DCSVAL, DNORIN, EV_fun, q_menu_fun, interpolate

!write(6, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') X, b_sup_short, b_inf_short, b_sup
b_short = MAX(MIN(X(1), b_sup_short), b_inf_short)

b_sup_local = b_sup !NO NEED TO CALCULATE ISSUANCE LIMIT
b_long  = MAX(MIN(X(2), b_sup_local), b_inf)


!PENALIZE IF DEBT IS OUTSIDE SPECIFIED RANGE. THIS PENALTY SPEEDS UP OPTIMIZATION ROUTINE. ALTERNATIVE TO HAVE
!A FLAT FUNCTION IN THAT REGION MAKES OPTIMIZATION ROUTINE EXAUST ITS MAXIMUM NUMBER OF ITERATIONS
penalty = 10d+0*(max(b_inf-X(2),0d+0)**2d+0 + max(X(2)-b_sup_local,0d+0)**2d+0 + max(b_inf_short-X(1),0d+0)**2d+0 + max(X(1)-b_sup_short,0d+0)**2d+0)

!CALCULATE EXPECTED CONTINUATION VALUE
exp_v_next = EV_fun(b_short, b_long, y_initial, i_access_initial)


d  = i_default_global
!q =  q_fun(b_short, b_long, y_initial)
q = q_menu_fun(b_short, b_long, y_initial, i_access_initial)
q_short = exp(-r)
borrowing =  b_long - b_long_initial*(1d+0-default_grid(d))*(1d+0 - delta)
output = EXP(y_initial)  

consum = output + b_short_initial - coupon*b_long_initial + borrowing*q - b_short*q_short - absortion

util = (max(consum, 1d-6))**(1-sigma)/(1-sigma)  !COMPUTE CURRENT UTILITY FLOW

objective_function =  - util - beta*exp_v_next + penalty

c_global = consum
util_global = util
ev_global = exp_v_next

end function


!RHS of BELLMAN EQUATION WHEN EXCLUDED FROM CREDIT MARKETS. NO DEBT IS CHOSEN OR REPAID AND THE ONLY CHOICE IS THE
!RESERVE ACCUMULATION
DOUBLE PRECISION function objective_function_excl(reserves)
USE param
DOUBLE PRECISION, INTENT(IN) :: reserves
INTEGER :: other_type, i,j,h, t, d
DOUBLE PRECISION :: q_fun, acum, exp_v_excl, exp_v_excl_ends, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, q, DNORDF, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, cdf, cdf1,&
                    cdf_threshold, DNORIN, q_risk_free, output, b_short, b_long, q_short, exp_v_next,&
                    EV_excl_fun, tax, tfp, labor, g, consum, util, interpolate, net_revenue, tax0_fun

EXTERNAL q_fun, DNORDF, y_fun, v0_fun, v1_fun, DCSVAL, DNORIN, EV_excl_fun, interpolate

b_long = b_long_initial*exp(r)
b_short = MAX(MIN(reserves, b_sup_short), b_inf_short)


q_short = exp(-r)
output = exp(y_initial)*(1-output_cost)
consum = output + b_short_initial - b_short*q_short - absortion

util = (max(consum, 1d-6))**(1-sigma)/(1-sigma)  !COMPUTE CURRENT UTILITY FLOW

exp_v_next = EV_excl_fun(b_short, b_long, y_initial, i_access_initial) !COMPUTE CONTINUATION VALUE UNDER DEFAULT

objective_function_excl = - util - beta * exp_v_next 


c_global = consum
output_global = output
end



!RHS of BELLMAN EQUATION FOR A GIVEN RESERVE LEVEL. THIS FUNCTION IS USED TO COMPUTE THE OPTIMAL DEBT LEVEL
!FOR A GIVEN CHOICE OF RESERVES
DOUBLE PRECISION function objective_function_debt(debt)
USE param
DOUBLE PRECISION, INTENT(IN) :: debt
INTEGER :: other_type, i,j,h, t, d
DOUBLE PRECISION :: q_fun, acum, exp_v_excl, exp_v_excl_ends, value_next, y_next, y_fun, scalar, DCSVAL, &
                    y_threshold, y_max, y_min, y_next1, acum2,acum1, q, DNORDF, v0_fun, v1_fun, &
                    value_next_tirar, acum_int, acum_int_lo, acum_int_hi, borrowing, cdf, cdf1,&
                    cdf_threshold, DNORIN, q_risk_free, output, b_short, b_long, q_short, exp_v_next, EV_fun,&
                    q_menu_fun, objective_function, penalty

EXTERNAL q_fun, DNORDF, y_fun, v0_fun, v1_fun, DCSVAL, DNORIN, EV_fun, q_menu_fun, objective_function

b_short = b_short_global_optimize
b_long  = MAX(MIN(debt, b_sup), b_inf)


objective_function_debt = objective_function((/b_short, b_long/))

end

    
subroutine optimize(b_next, v_value)
USE param
integer :: t, d, index_opt, index_vector(1), num, num_tirar, N, IPARAM(7), IBTYPE,IERSVR, IPACT, ISACT,&
           index_short, index_long, i_access, i, ITER, NP, num_long
parameter (num_long=100, num_tirar = 200, N=2, NP = 2)
DOUBLE PRECISION :: b_next(2), v_value, new_value, q_fun, old_value, b, q,vector(b_num_long),&
                    b_next_short, b_next_long, b_next_guess(2), b_tirar(2), b_next_inf, b_next_sup, &
                    difference, XSCALE(N), FSCALE, RPARAM(7), FVALUE, XUB(N), XLB(N), fun_right, fun_left,&
                    derivative, eps,  y_threshold, y_fun, y_current_type_1, y_current_type_2,DNORDF, q_menu_fun,&
                    objective_function, b_short_next_inf, b_short_next_sup, FTOL, matrix_direction(N,N), &
                    b_next_graph(2)

EXTERNAL q_fun, DUVMIF, objective_function, ERSET, y_fun, DNORDF, q_menu_fun


b_next_inf = b_inf
b_next_sup = b_sup

if (indicator_global_search > 0) then  !SEARCH FOR BEST INITIAL GUESS USING ALL GRID RANGE

   old_value = 10d+15
   index_opt = 1
   do i=1,num_short_next

      b_next_short = b_next_short_grid(i)
      b_short_global_optimize = b_next_short
   !FIND THE OPTIMAL DEBT ACCUMULATION AND IMPLIED CONTINUATION VALUE FOR A CHOICE OF RESERVES = b_next_short
      call optimize_debt(b_next_long, new_value)
       if (new_value <= old_value) then
          index_opt = i
           b_next_guess(1) = b_next_short
           b_next_guess(2) = b_next_long
           old_value = new_value
           index_short = i
       end if
   end do

else !FOCUS GRID SEARCH IN A NEIGHBORHOOD OF PREVIOUS OPTIMAL VALUE

b_next_guess(1) = b0_next_short_matrix(i_b_short_global, i_b_long_global, i_y_global, i_access_initial)
b_next_guess(2) = b0_next_long_matrix(i_b_short_global, i_b_long_global, i_y_global, i_access_initial)
index_short = 2 !SET IT TO A NUMBER > 1 & < num_short_next SO THAT THE OPTIMIZATION ROUTINE IS INVOKED


end if

XSCALE = 1d+0
FSCALE = 1d+0
IPARAM(1) = 0d+0
IBTYPE = 0
XLB(1) = b_inf_short
XLB(2) = b_inf
XUB(1) = b_sup_short
XUB(2) = b_sup


if (index_short ==1) then
    b_next = b_next_guess
    v_value = -objective_function(b_next)
elseif(index_short ==num_short_next ) then
    b_next = b_next_guess
    v_value = -objective_function(b_next)
else


    matrix_direction = 0d+0
    matrix_direction(1,1) = 1d+0
    matrix_direction(2,2) = 1d+0
    FTOL = 1d-8

    call POWELL(objective_function, b_next_guess, matrix_direction, N, NP, FTOL, ITER, FVALUE)
    b_next(1) = MAX(MIN(b_next_guess(1), b_sup_short), b_inf_short)
    b_next(2)  = MAX(MIN(b_next_guess(2), b_sup), b_inf)
    v_value = -FVALUE

end if



10 end subroutine


subroutine optimize_excl(b_next, v_value)
USE param
integer :: MAXFN, t, d, index_opt, index_vector(1), num, num_tirar, ITER, NP, N, i
parameter (num=20, num_tirar = 500, N = 1, NP =1)
DOUBLE PRECISION :: b_next, v_value, objective_function_excl, new_value, q_fun, old_value, b, q,&
                    b_next_grid(num), STEP, BOUND, XACC, b_next_guess, b_tirar, b_next_inf, b_next_sup,&
                    b_tirar_inf, b_tirar_sup
EXTERNAL objective_function_excl, DUVMIF
b_next_inf = b_inf_short
b_next_sup = b_sup_short


old_value = 10d+3
do i=1,num
   b_next_grid(i) = b_next_inf + (b_next_sup - b_next_inf) * (i-1)/ (num - 1)
   new_value = objective_function_excl(b_next_grid(i))
   if (new_value <= old_value) then
      index_opt = i
      b_next_guess = b_next_grid(i)
      old_value = new_value
  end if
end do


STEP = (b_next_grid(2) - b_next_grid(1))*0.5
BOUND = 10*EXP(y_sup)
XACC = 1d-6
MAXFN = 1000


b_tirar_inf = b_inf_short + 1.0d-4
b_tirar_sup = b_sup_short - 1.0d-4
if (index_opt == 1 .AND. objective_function_excl(b_tirar_inf) > objective_function_excl(b_next_grid(index_opt))) then
   b_next = b_next_guess
   v_value = -old_value
ELSEIF(index_opt == 1 .AND. objective_function_excl(b_tirar_inf) < objective_function_excl(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO 0 BUT NOT AT ZERO
       b_next_guess = b_tirar_inf
       call DUVMIF(objective_function_excl, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = -objective_function_excl(b_next)

elseif (index_opt == num .AND. objective_function_excl(b_tirar_sup) >= objective_function_excl(b_next_grid(index_opt))) then
   b_next = b_next_guess
   v_value = -old_value
elseif(index_opt == num .AND. objective_function_excl(b_tirar_sup) < objective_function_excl(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MAXIMUM GRID POINT BUT NOT AT MAXIMUM GRID POINT
       b_next_guess = b_tirar_sup
       call DUVMIF(objective_function_excl, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = -objective_function_excl(b_next)

else
   if (ABS(old_value - objective_function_excl(b_next_grid(index_opt-1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10 .OR.  &
       ABS(old_value - objective_function_excl(b_next_grid(index_opt+1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10) then
!function is too flat to invoke optimization routine
       b_next = b_next_guess
       v_value = -old_value
   else
      call DUVMIF(objective_function_excl, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = -objective_function_excl(b_next)
   end if
end if
 


10 end subroutine


subroutine optimize_debt(b_next, v_value)
USE param
integer :: MAXFN, t, d, index_opt, index_vector(1), num, i
parameter (num=15)
DOUBLE PRECISION :: b_next, v_value, objective_function_debt, new_value, q_fun, old_value, b, q, q_menu_fun,&
                    b_next_grid(num), STEP, BOUND, XACC, b_next_guess, b_tirar_sup, b_tirar_inf, b_next_inf, b_next_sup,&
                    ancho, debt_inf, debt_sup
EXTERNAL objective_function_debt, DUVMIF, q_menu_fun

!open(100, FILE ='tirar.txt',STATUS='replace')

b_next_inf = b_inf
b_next_sup = b_sup

if (indicator_global_search > 0) then  !SEARCH FOR BEST INITIAL GUESS USING ALL GRID RANGE
  old_value = 10d+15
  index_opt = 1
   do i=1,num
    !DEFINE THE GRID IN A NEIGHBORHOOD OF THE PREVIOUS OPTIMAL VALUE
!    b_next_grid(i) = debt_inf + (debt_sup - debt_inf) * (i-1)/ (num - 1) !DEFINE GRID BASED ON PREVIOUS OPTIMAL VALUE
     b_next_grid(i) = b_next_inf + (b_next_sup - b_next_inf) * (i-1)/ (num - 1)
     new_value = objective_function_debt(b_next_grid(i))
     if (new_value <= old_value) then
        index_opt = i
        b_next_guess = b_next_grid(i)
        old_value = new_value
     end if
   end do
else !FOCUS GRID SEARCH IN A NEIGHBORHOOD OF PREVIOUS OPTIMAL VALUE


   ancho = 0.25d+0
!IF POLICY FUNCTIONS HAVE CONVERGED, USE PREVIOUS OPTIMAL VALUE FOR DEBT TO NARROW GRID SEARCH
 debt_inf = MAX(b_next_inf, b0_next_long_matrix(i_b_short_global, i_b_long_global, i_y_global, i_access_initial) - ancho)
 debt_sup = MIN(b_next_sup, b0_next_long_matrix(i_b_short_global, i_b_long_global, i_y_global, i_access_initial) + ancho)

  old_value = 10d+15
  index_opt = 1
  do i=1,num
    !DEFINE THE GRID IN A NEIGHBORHOOD OF THE PREVIOUS OPTIMAL VALUE
    b_next_grid(i) = debt_inf + (debt_sup - debt_inf) * (i-1)/ (num - 1) !DEFINE GRID BASED ON PREVIOUS OPTIMAL VALUE
    new_value = objective_function_debt(b_next_grid(i))
    if (new_value <= old_value) then
       index_opt = i
       b_next_guess = b_next_grid(i)
       old_value = new_value
    end if
!write(nout, '(I3, X, F12.8, X, I3, X, F12.8, X, F12.8)') i, b_next_grid(i), index_opt, new_value, &
!q_menu_fun(b_short_global_optimize, b_next_grid(i), y_initial, i_access_initial)

!write(100, '(F12.8, X, F12.8, X, F12.8)') b_next_grid(i), MIN(new_value, 100)
  end do
end if

!CLOSE(100)
!pause
STEP = (b_next_grid(2) - b_next_grid(1))*0.5
BOUND = 10*EXP(y_sup)
XACC = 1d-5
MAXFN = 1000



b_tirar_sup = b_next_sup - 1.0d-4
b_tirar_inf = b_next_inf + 1.0d-4

if (index_opt == 1 .AND. objective_function_debt(b_tirar_inf) >= objective_function_debt(b_next_grid(index_opt))) then
   b_next = b_next_guess
   v_value = old_value

elseif(index_opt == 1 .AND. objective_function_debt(b_tirar_inf) < objective_function_debt(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MINIMUM GRID POINT BUT NOT AT MINIMUM GRID POINT
       b_next_guess = b_tirar_inf
       call DUVMIF(objective_function_debt, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = objective_function_debt(b_next)
elseif (index_opt == num .AND. objective_function_debt(b_tirar_sup) >= objective_function_debt(b_next_grid(index_opt))) then
       b_next = b_next_guess
       v_value = old_value
elseif(index_opt == num .AND. objective_function_debt(b_tirar_sup) < objective_function_debt(b_next_grid(index_opt))) then
       !THE MAXIMUM MAY BE CLOSE TO MAXIMUM GRID POINT BUT NOT AT MAXIMUM GRID POINT
       b_next_guess = b_tirar_sup
       call DUVMIF(objective_function_debt, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = objective_function_debt(b_next)
else
   if (ABS(old_value - objective_function_debt(b_next_grid(index_opt-1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10 .OR.  &
       ABS(old_value - objective_function_debt(b_next_grid(index_opt+1)))/(b_next_grid(2)-b_next_grid(1)) < 1d-10) then
       b_next = b_next_guess
       v_value = old_value
   else
      call DUVMIF(objective_function_debt, b_next_guess, STEP, BOUND, XACC, MAXFN, b_next)
       v_value = objective_function_debt(b_next)
   end if
end if

10 end subroutine

DOUBLE PRECISION function interpolate(net, y, matrix)
USE param
DOUBLE PRECISION, INTENT(IN) :: net, y, matrix(net_num, y_num_finer)
DOUBLE PRECISION ::  slope, weight, acum, ratio_net, ratio_y
INTEGER :: index_net, index_y, i_net, i_y



index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))
index_net = (MAX(MIN(INT((net_num-1)*(net-net_grid(1))/(net_grid(net_num)-net_grid(1)))+1,net_num-1), 1))


ratio_net = (net - net_grid(index_net)) / (net_grid(index_net+1) - net_grid(index_net))
ratio_y = (y - y_grid_finer(index_y)) / (y_grid_finer(index_y+1) - y_grid_finer(index_y))
acum=0



do i_net=0,1
   do i_y=0,1
        weight = ((1-i_net)*(1 - ratio_net) + i_net*ratio_net )* ((1-i_y)*(1 - ratio_y) + i_y*ratio_y )
        acum = acum + matrix(index_net + i_net, index_y + i_y)*weight
!        if (indicator_tirar>0) then
!            WRITE(nout, '(I4, X, I4, X, F12.8, X, F12.8, X, F12.8, X, F12.8, X, I4, X, F12.8)') &
!            index_net + i_net, index_y + i_y, net_grid(index_net + i_net), y_grid(index_y + i_y), matrix(index_net + i_net, index_y + i_y), acum
!            i_y*ratio_y, i_y,ratio_y
!        end if
   end do
end do
interpolate = acum
    end

    
!FUNCTION THAT APPROXIMATES THE blong LEVEL AT WHICH q(bshort, blong, y, m) = qbar IF THERE IS SUCH VALUE
!FOR blong IN [b_inf, b_sup]
DOUBLE PRECISION function b_limit_fun_spline(b_short, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DCSVAL, value_left, value_right
INTEGER :: index_b, NINTV, index_y
EXTERNAL DCSVAL

NINTV = b_num_short_finer - 1
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))

value_left  = DCSVAL(b_short, NINTV, break_matrix_b_limit(:,index_y, index_access), coeff_matrix_b_limit(:,:,index_y, index_access))
value_right = DCSVAL(b_short, NINTV, break_matrix_b_limit(:,index_y+1, index_access), coeff_matrix_b_limit(:,:,index_y+1, index_access))


slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
b_limit_fun_spline = value_left + slope*(y - y_grid_finer(index_y))

end
    
    
    
    DOUBLE PRECISION function v0_fun(b_short, b_long, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL

!NINTV = b_num - 1
index_y = (MAX(MIN(INT((y_num-1)*(y-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))

!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_short+1d-6, b_short),b_sup_short-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

value_left  = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots, b_long_knots, b_num_short, b_num_long, &
coeff_matrix_v0(:,:,index_y, index_access))

value_right = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots, b_long_knots, b_num_short, b_num_long, &
coeff_matrix_v0(:,:,index_y+1, index_access))

slope = (value_right - value_left)/ (y_grid(index_y+1) - y_grid(index_y))

v0_fun = value_left + slope * (y - y_grid(index_y))

    end function

   

DOUBLE PRECISION function v1_fun(b_short, b_long, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_b, NINTV
EXTERNAL DBS2VL, interpolate3

v1_fun = interpolate3(b_short, b_long, y, v1_matrix(:,:,:, index_access))

end function

    


DOUBLE PRECISION function interpolate3(bs, bl, y, matrix)
USE param
DOUBLE PRECISION, INTENT(IN) :: bs, bl, y, matrix(b_num_short, b_num_long, y_num)
DOUBLE PRECISION ::  slope, weight, acum, ratio_bs, ratio_bl, ratio_y, eps
INTEGER :: index_bs, index_bl, index_y, i_bs, i_bl, i_y

eps = 1d-8 !NEEDED TO ADD EPSILON BECAUSE IF y = y_grid(i) IT MAY CHOOSE index_y = i-1 DUE TO APPROX ERROR IN THE ROUNDING
index_y = (MAX(MIN(INT((y_num-1)*(y+eps-y_grid(1))/(y_grid(y_num)-y_grid(1)))+1,y_num-1), 1))
index_bs = (MAX(MIN(INT((b_num_short-1)*(bs+eps-b_grid_short(1))/(b_grid_short(b_num_short)-b_grid_short(1)))+1,b_num_short-1), 1))
index_bl = (MAX(MIN(INT((b_num_long -1)*(bl+eps-b_grid_long(1)) /(b_grid_long(b_num_long)-b_grid_long(1)))+1,b_num_long-1), 1))



ratio_bs = (bs - b_grid_short(index_bs)) / (b_grid_short(index_bs+1) - b_grid_short(index_bs))
ratio_bl = (bl - b_grid_long(index_bl)) / (b_grid_long(index_bl+1) - b_grid_long(index_bl))
ratio_y = (y - y_grid(index_y)) / (y_grid(index_y+1) - y_grid(index_y))
acum=0

do i_bs=0,1
   do i_bl = 0,1
      do i_y=0,1
        weight = ((1-i_bs)*(1 - ratio_bs) + i_bs*ratio_bs) * ((1-i_bl)*(1 - ratio_bl) + i_bl*ratio_bl) * &
                 ((1-i_y)*(1 - ratio_y) + i_y*ratio_y )
        acum = acum + matrix(index_bs + i_bs, index_bl + i_bl, index_y + i_y)*weight
      end do
   end do
end do

interpolate3 = acum
end function


DOUBLE PRECISION function q_paid_fun(b_short, b_long, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl, interpolate3
INTEGER :: index_b, NINTV
EXTERNAL DBS2VL, interpolate3

q_paid_fun = interpolate3(b_short, b_long, y, q_nodef_matrix(:,:,:, index_access))

end


    
    
    
    
DOUBLE PRECISION function EV_fun(b_short, b_long, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL

!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_short+1d-6, b_short),b_sup_short-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

!NINTV = b_num - 1
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))

value_left  = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots_finer, b_long_knots_finer, b_num_short_finer, b_num_long_finer, &
coeff_matrix_ev(:,:,index_y, index_access))

value_right = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots_finer, b_long_knots_finer, b_num_short_finer, b_num_long_finer, &
coeff_matrix_ev(:,:,index_y+1, index_access))

slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
EV_fun= value_left + slope * (y - y_grid_finer(index_y))

end function


    
    
DOUBLE PRECISION function EV_excl_fun(b_short, b_long, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL

!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_short+1d-6, b_short),b_sup_short-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

!NINTV = b_num - 1
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))

value_left  = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots_finer, b_long_knots_finer, b_num_short_finer, b_num_long_finer, &
coeff_matrix_ev_excl(:,:,index_y, index_access))

value_right = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots_finer, b_long_knots_finer, b_num_short_finer, b_num_long_finer, &
coeff_matrix_ev_excl(:,:,index_y+1, index_access))

slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
EV_excl_fun= value_left + slope * (y - y_grid_finer(index_y))

end function

    
    

DOUBLE PRECISION function q_menu_fun(b_short, b_long, y, index_access)
USE param
DOUBLE PRECISION, INTENT(IN) :: b_short, b_long, y
INTEGER, INTENT(IN) :: index_access
DOUBLE PRECISION ::  slope, DBS2VL, value_left, value_right, bs, bl
INTEGER :: index_b, NINTV, index_y
EXTERNAL DBS2VL

!BOUNDARY VALUES ARE ADJUSTED SO THAT bs AND bl ARE ALWAYS WITHING THE GRID REGION
bs=  min(max(b_inf_short+1d-6, b_short),b_sup_short-1d-6)
bl=  min(max(b_inf+1d-6, b_long),b_sup-1d-6)

!NINTV = b_num - 1
index_y = (MAX(MIN(INT((y_num_finer-1)*(y-y_grid_finer(1))/(y_grid_finer(y_num_finer)-y_grid_finer(1)))+1,y_num_finer-1), 1))


value_left  = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots_finer, b_long_knots_finer, b_num_short_finer, b_num_long_finer, &
coeff_matrix_q_menu(:,:,index_y, index_access))

value_right = DBS2VL(bs,bl, KORDER, KORDER, b_short_knots_finer, b_long_knots_finer, b_num_short_finer, b_num_long_finer, &
coeff_matrix_q_menu(:,:,index_y+1, index_access))

slope = (value_right - value_left)/ (y_grid_finer(index_y+1) - y_grid_finer(index_y))
q_menu_fun  = value_left + slope * (y - y_grid_finer(index_y))

end


subroutine compute_knots
USE param
!INTEGER :: KORDER
!PARAMETER(KORDER=3)
!DOUBLE PRECISION, DIMENSION (b_num_long + KORDER) :: b_long_knots
!DOUBLE PRECISION, DIMENSION (b_num_short + KORDER) :: b_short_knots
EXTERNAL DBSNAK

CALL DBSNAK (b_num_long, b_grid_long, 3, b_long_knots)
CALL DBSNAK (b_num_short, b_grid_short,3, b_short_knots)

CALL DBSNAK (b_num_long_finer, b_grid_long_finer, 3, b_long_knots_finer)
CALL DBSNAK (b_num_short_finer, b_grid_short_finer,3, b_short_knots_finer)

end subroutine


double precision function b_limit_fun(bnext)
use param
double precision, intent(in) :: bnext
double precision :: q_menu_fun
external q_menu_fun

b_limit_fun =  q_menu_fun(b_short_next_global, bnext, y_initial, i_access_initial) - q_global
end function  
    
    
!COMPUTE MATRICES AND SPLINE COEFFICIENTS FOR q(bs, bl, y) and E[V(bs, bl, y') | y)]
subroutine compute_q_ev
USE param
INTEGER :: i_b_short, i_b_long, i_y, i, ILEFT, IRIGHT, i_access, index_access_next, KXORD, KYORD, LDF, MAXFN
parameter (KXORD =3, KYORD = 3)
DOUBLE PRECISION :: b_short, b_long, y_threshold, cdf_threshold, y_fun, DNORDF, exp_v_next, scalar, y_min, y_max, DNORIN, &
                    acum1, acum2, acum(2), cdf, cdf1, y_next, y_next1, value_next, v0_fun, v1_fun, q_fun, exp_v_excl, &
                    exp_v_excl_ends, FDATA(b_num_short_finer), DLEFT, DRIGHT, BREAK_GRID(b_num_short_finer), &
                    CSCOEF(4,b_num_short_finer), DCSVAL, tirar, FDATA_2D_finer(b_num_short_finer, b_num_long_finer), &
                    BSCOEF_finer(b_num_short_finer, b_num_long_finer), q_menu_fun,&
                    left, right, ERRREL, ERRABS, b_limit_fun, b_limit_fun_spline

DOUBLE PRECISION, DIMENSION (b_num_short_finer, b_num_long_finer, y_num_finer, 2) :: EV_matrix, q_menu_matrix, EV_excl_matrix
DOUBLE PRECISION, DIMENSION (b_num_short_finer, y_num_finer, 2) :: blong_limit_matrix

EXTERNAL y_fun, DNORDF, v0_fun, v1_fun, DCSDEC, DNORIN, q_fun, DBS2IN, q_menu_fun, DZBREN, b_limit_fun, b_limit_fun_spline


do i_access = 1,2
  i_access_initial = i_access
  do i_y = 1,y_num_finer
 ! WRITE(nout, *) i_y
    y_initial = y_grid_finer(i_y)
      do i_b_short = 1, b_num_short_finer
      b_short = b_grid_short_finer(i_b_short)
     
          do i_b_long = 1, b_num_long_finer
          b_long = b_grid_long_finer(i_b_long)

            do index_access_next = 1,access_num

 !COMPUTE EXPECTED VALUE FUNCTION FOR THE NEXT PERIOD IF IT CHOOSES A PORTFOLIO (b_short, b_long)
 !WHEN CURRENT INCOME = y_initial
               y_threshold   = y_fun(b_short, b_long, index_access_next)
               cdf_threshold = DNORDF((y_threshold - rho*y_initial - (1-rho)* mean_y)/std_eps)

               exp_v_next = 0d+0
               scalar = 1d+0 / SQRT(pi_number)

	       y_min = (1-rho)*mean_y + rho* y_initial - 4*std_eps
               y_max = (1-rho)*mean_y + rho* y_initial + 4*std_eps

	       acum1=0d+0
	       acum2=0d+0


              if (cdf_threshold <= cdf_inf) then
                 do i=1,quad_num_v
          	    cdf = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
          	    y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
          	    value_next = v0_fun(b_short, b_long, y_next, index_access_next)
          	    acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
       	         end do
              elseif(cdf_threshold >= cdf_sup) then

       	         do i=1,eps_num
          	    cdf = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
          	    y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
          	    value_next = v1_fun(b_short, b_long, y_next, index_access_next)
             	    acum1 = acum1 + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
       	         end do

   	      else
      	         do i=1,eps_num

        	    cdf = 0.5*(quad_x_v(i) + 1)*(cdf_threshold - cdf_inf) + cdf_inf
        	    y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
         	    value_next = v1_fun(b_short, b_long, y_next, index_access_next)
        	    acum1 = acum1 + 0.5*(cdf_threshold - cdf_inf)*quad_w_v(i) * value_next

	            cdf1 = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_threshold) + cdf_threshold
        	    y_next1 = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf1) * std_eps !DCSVAL(cdf1, NINTV, BREAK_eps, CSCOEF_eps)
	            value_next = v0_fun(b_short, b_long, y_next1, index_access_next)
               	    acum2 = acum2 + 0.5*(cdf_sup - cdf_threshold)*quad_w_v(i) * value_next
!           WRITE(119, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') y_next, v1_fun(y_next), y_next1, value_next
      	         end do
    	      end if
            acum(index_access_next) = trans_matrix(i_access, index_access_next) * (acum1+acum2)/(cdf_sup - cdf_inf)
            end do

           EV_matrix(i_b_short, i_b_long, i_y, i_access) = SUM(acum)
           q_menu_matrix(i_b_short, i_b_long, i_y, i_access) = q_fun(b_short, b_long, y_initial, i_access)
        
          
          
 !COMPUTE EXPECTED VALUE FUNCTION FOR THE NEXT PERIOD IF EXCLUDED AND CHOOSES A PORTFOLIO (b_short)
 !WHEN CURRENT INCOME = y_initial

          acum = 0d+0
       	   do index_access_next = 1,access_num
    	      exp_v_excl      = 0d+0
	          exp_v_excl_ends = 0d+0
	          scalar = 1d+0 / SQRT(pi_number)

               do i=1,quad_num_v
                  !CONTINUES EXCLUDED IN THE NEXT PERIOD (NO BORROWING AND DEFAULT COST)
                  cdf = 0.5*(quad_x_v(i) + 1)*(cdf_sup - cdf_inf) + cdf_inf
                  y_next  = (1-rho)*mean_y + rho* y_initial + DNORIN(cdf) * std_eps
                  value_next = v1_fun(b_short, b_long, y_next, index_access_next)
                  exp_v_excl = exp_v_excl + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next

                  !EXCLUSION ENDS IN THE NEXT PERIOD (STARTS WITH RECOVERY DEBT b_recov AND DECIDES WHETHER TO DEFAULT)
                  value_next = max(v0_fun(b_short, min(b_recov, b_long), y_next, index_access_next), v1_fun(b_short, min(b_recov, b_long), y_next, index_access_next))
                  exp_v_excl_ends = exp_v_excl_ends + 0.5*(cdf_sup - cdf_inf)*quad_w_v(i) * value_next
                  
                  !if(i_access ==2 .and. i_b_short ==1 .and. i_b_long ==1) then
                  !    write(6, '(F12.8, X, F14.8, X, F14.8)') y_next, v1_fun(b_short, b_long, y_next, index_access_next), max(v0_fun(b_short, b_long, y_next, index_access_next), v1_fun(b_short, b_long, y_next, index_access_next))
                  !end if
                  

               end do
               acum(index_access_next) = trans_matrix(i_access, index_access_next) * &
               ((1d+0 - prob_excl_end) * exp_v_excl + prob_excl_end * exp_v_excl_ends)/(cdf_sup - cdf_inf)
             end do
          
                 EV_excl_matrix(i_b_short, i_b_long, i_y, i_access) = SUM(acum)
                 !if(i_access ==2 .and. i_b_short ==1 .and. i_b_long ==1) then
                 !    write(nout, '(I4, X, F12.8, X, F12.8)') i_y, y_initial, EV_excl_matrix(i_b_short, i_b_long, i_y, i_access) 
                 !    pause       
                 !end if
                 !
          
          end do

 !          WRITE(nout, *) i_y, i_b_short, EV_excl_matrix(i_b_short, i_y, i_access)
      end do
  end do
end do



!COMPUTE COEFFICIENTS OF MATRIX 
do i_access = 1,access_num
   do i_y = 1,y_num_finer

       FDATA_2D_finer = EV_matrix(:,:,i_y,i_access) 
       LDF = size(FDATA_2D_finer,1)
    !write(6,*) 'cesar', cesar
       CALL DBS2IN (b_num_short_finer, b_grid_short_finer, b_num_long_finer, b_grid_long_finer, FDATA_2D_finer, LDF, KXORD, KYORD, &
                    b_short_knots_finer, b_long_knots_finer, BSCOEF_finer)
       coeff_matrix_ev(:,:,i_y,i_access) = BSCOEF_finer

       
       FDATA_2D_finer = q_menu_matrix(:,:,i_y,i_access) 
       LDF = size(FDATA_2D_finer,1)
    !write(6,*) 'cesar', cesar
       CALL DBS2IN (b_num_short_finer, b_grid_short_finer, b_num_long_finer, b_grid_long_finer, FDATA_2D_finer, LDF, KXORD, KYORD, &
                    b_short_knots_finer, b_long_knots_finer, BSCOEF_finer)
       coeff_matrix_q_menu(:,:,i_y,i_access) = BSCOEF_finer


       FDATA_2D_finer = EV_excl_matrix(:,:,i_y,i_access) 
       LDF = size(FDATA_2D_finer,1)
    !write(6,*) 'cesar', cesar
       CALL DBS2IN (b_num_short_finer, b_grid_short_finer, b_num_long_finer, b_grid_long_finer, FDATA_2D_finer, LDF, KXORD, KYORD, &
                    b_short_knots_finer, b_long_knots_finer, BSCOEF_finer)
       coeff_matrix_ev_excl(:,:,i_y,i_access) = BSCOEF_finer

       

   end do
end do




!USE THE SPLINE COEFFICIENTS FOR q_menu_fun TO CALCULATE THE blong LEVELS AT WHICH q(bshort, blong, y, m) = qbar IF THERE IS SUCH VALUE
!FOR blong IN [b_inf, b_sup]


do i_access = 1,access_num
  i_access_initial = i_access
  do i_y = 1,y_num_finer
    y_initial = y_grid_finer(i_y)
      do i_b_short = 1, b_num_short_finer
        b_short = b_grid_short_finer(i_b_short)
        b_short_next_global = b_short 
        if (q_menu_fun(b_short, b_sup, y_initial, i_access_initial) >= q_global) then 
            blong_limit_matrix(i_b_short, i_y, i_access) = b_sup
        elseif  (q_menu_fun(b_short, b_sup, y_initial, i_access_initial) < q_global .and. &
            q_menu_fun(b_short, b_inf, y_initial, i_access_initial) > q_global) then
            !MAXIMUM DEBT > b*(1-delta). CAN ISSUE LONG-TERM DEBT UNTIL q(b)=qbar
            left = b_inf
            right = b_sup
            MAXFN=1000
            ERRREL = 1D-6
            ERRABS = 1D-6
            CALL DZBREN (b_limit_fun, ERRABS, ERRREL, left, right, MAXFN)
            blong_limit_matrix(i_b_short, i_y, i_access) = right
        else
            !MAXIMUM DEBT = b_inf. CAN BUY BACK DEBT BUT CANNOT ISSUE NEW LONG-TERM DEBT
            blong_limit_matrix(i_b_short, i_y, i_access) = b_inf
        end if

        !if(i_access == 2 ) then
        !    write(6, '(F12.8, X, F12.8, X, F12.8, X, F12.8)') q_menu_fun(b_short, b_inf, y_initial, i_access_initial), q_menu_fun(b_short, b_sup, y_initial, i_access_initial), &
        !    blong_limit_matrix(i_b_short, i_y, i_access), q_menu_fun(b_short, blong_limit_matrix(i_b_short, i_y, i_access), y_initial, i_access_initial)
        !    pause
        !end if
        
      end do
  end do
end do


do i_access = 1,access_num
   do i_y = 1,y_num_finer
        FDATA = blong_limit_matrix(:,i_y,i_access) 
	    ILEFT  = 0
	    IRIGHT = 0
	    call  DCSDEC (b_num_short_finer, b_grid_short_finer, FDATA, ILEFT, DLEFT, IRIGHT, DRIGHT, BREAK_GRID, CSCOEF)
        break_matrix_b_limit(:,i_y, i_access) = BREAK_GRID
        coeff_matrix_b_limit(:,:,i_y, i_access) = CSCOEF
      
   end do
end do


open(100, FILE ='coeff_excl.txt',STATUS='replace')
open(101, FILE ='coeff_EV.txt',STATUS='replace')
open(102, FILE ='coeff_q.txt',STATUS='replace')
open(103, FILE ='coeff_blimit.txt',STATUS='replace')
open(104, FILE ='EV.txt',STATUS='replace')

do i_access = 1,access_num
   do i_y = 1,y_num_finer
       do i_b_short = 1,b_num_short_finer
          write(103, '(F20.8, X, F20.8, X, F20.8, X, F20.8, X, F20.8)') break_matrix_b_limit(i_b_short,i_y, i_access), coeff_matrix_b_limit(1,i_b_short,i_y, i_access), &
          coeff_matrix_b_limit(2,i_b_short,i_y, i_access), coeff_matrix_b_limit(3,i_b_short,i_y, i_access), coeff_matrix_b_limit(4,i_b_short,i_y, i_access)
          do i_b_long = 1,b_num_long_finer
           
             WRITE(100, '(F20.8)') coeff_matrix_EV_excl(i_b_short, i_b_long, i_y, i_access)
             WRITE(101, '(F20.8)') coeff_matrix_EV(i_b_short, i_b_long, i_y, i_access)
             write(102, '(F20.8)') coeff_matrix_q_menu(i_b_short, i_b_long, i_y, i_access)
             write(104, '(F20.8)') EV_matrix(i_b_short, i_b_long, i_y, i_access)
          end do
       end do
   end do
end do
CLOSE(100)
CLOSE(101)
CLOSE(102)
CLOSE(103)
close(104)




    end subroutine



subroutine read_q_ev
USE param
INTEGER :: i_b_short, i_b_long, i_y, i, i_access, index_access_next

open(100, FILE ='coeff_excl.txt')
open(101, FILE ='coeff_EV.txt')
open(102, FILE ='coeff_q.txt')
open(103, FILE ='coeff_blimit.txt')
do i_access = 1,access_num
   do i_y = 1,y_num_finer
       do i_b_short = 1,b_num_short_finer
          read(103, '(F20.8, X, F20.8, X, F20.8, X, F20.8, X, F20.8)') break_matrix_b_limit(i_b_short,i_y, i_access), coeff_matrix_b_limit(1,i_b_short,i_y, i_access), &
          coeff_matrix_b_limit(2,i_b_short,i_y, i_access), coeff_matrix_b_limit(3,i_b_short,i_y, i_access), coeff_matrix_b_limit(4,i_b_short,i_y, i_access)
          do i_b_long = 1,b_num_long_finer
             read(100, '(F20.8)') coeff_matrix_EV_excl(i_b_short, i_b_long, i_y, i_access)
             read(101, '(F20.8)')     coeff_matrix_EV(i_b_short, i_b_long, i_y, i_access)
             read(102, '(F20.8)') coeff_matrix_q_menu(i_b_short, i_b_long, i_y, i_access)
           end do
       end do
   end do
end do
CLOSE(100)
CLOSE(101)
CLOSE(102)
close(103)

end subroutine



subroutine iterate
USE param
DOUBLE PRECISION :: y_valor, b_valor, b_valor_def, v_valor, convergence, criteria, deviation, q_fun, b_next_short,&
                    b_next_long, g, b0_next(2), b1_next(2), b, q_value, v_valor_exl, FVALUE, q_menu_fun,&
                    v0_value, v1_value, q, acum_v, acum_excl, dev_q, dev_debt, dev_res, objective_function_excl, &
                    FDATA(b_num_short), DLEFT, DRIGHT, BREAK_GRID(b_num_short), CSCOEF(4,b_num_short), b_tirar, DCSVAL,&
                    FDATA_2D(b_num_short,b_num_long), BSCOEF(b_num_short, b_num_long), foc_vector(2)

DOUBLE PRECISION, DIMENSION(b_num_short, b_num_long,  y_num, 2) :: v0_matrix_new, v1_matrix_new, q_nodef_matrix_new, &
                                                                   tax0_matrix_new, tax1_matrix_new, output0_matrix, output1_matrix, &
                                                                   c0_matrix, c1_matrix, g0_matrix_new, g1_matrix_new

DOUBLE PRECISION, DIMENSION(b_num_short, b_num_long,  y_num, 2) :: v_matrix_new, dev_matrix, q_paid_matrix_new, dev_matrix_q, &
                             b0_next_short_matrix_new, b0_next_long_matrix_new, b1_next_short_matrix_new
INTEGER :: d, i, i_b_short, i_b_long, i_y, i_access, i_def_opt, i_b_zero, i_b_optimal_short0, i_b_optimal_long0, &
           i_b_optimal_short1, i_b_optimal_long1, i_b_next_short, i_b_next_long, ILEFT, IRIGHT,NINTV, N, &
           KXORD, KYORD, LDF, counter_final
parameter (KXORD =3, KYORD = 3)
EXTERNAL q_fun, DCSDEC, DCSVAL, objective_function_excl, q_menu_fun, DBS2IN



i_b_zero = b_num_long
criteria = 1d-6   !CRITERIA FOR CONVERGENCE
convergence = -1
indicator_global_search = 1
counter_final = 0

do WHILE(convergence<0)
   deviation = 0
   dev_q =0
   dev_res = 0
   dev_debt = 0
   counter = counter + 1
   counter_final = counter_final + 1
   q_global =  kappa* coupon*(1d+0-((1d+0-delta)*exp(-r))**counter)/(r+delta)

   do i_access = 1,access_num
      do i_y = 1,y_num

         FDATA_2D = v0_matrix(:,:,i_y,i_access) 
         LDF = size(FDATA_2D,1)
         CALL DBS2IN (b_num_short, b_grid_short, b_num_long, b_grid_long, FDATA_2D, LDF, KXORD, KYORD, b_short_knots, b_long_knots, BSCOEF)
         coeff_matrix_v0(:,:,i_y,i_access) = BSCOEF

         FDATA_2D = q_nodef_matrix(:,:,i_y,i_access) 
         LDF = size(FDATA_2D,1)
         CALL DBS2IN (b_num_short, b_grid_short, b_num_long, b_grid_long, FDATA_2D, LDF, KXORD, KYORD, b_short_knots, b_long_knots, BSCOEF)
         coeff_matrix_q(:,:,i_y,i_access) = BSCOEF

      end do
   end do

   call compute_q_ev !COMPUTE MATRICES OF SPLINE COEFFICIENTS TO APPROXIMATE EXPECTED VALUE FUNCTIONS AND BOND PRICE FUNCTION

      do i_access = 1,2
         i_access_initial = i_access
         g_initial = g_grid(i_access)
         do i_y = 1,y_num
            i_y_global = i_y
            y_initial = y_grid(i_y)

               do i_b_short =  1, b_num_short
                   i_b_short_global = i_b_short
                   b_short_initial = b_grid_short(i_b_short)
   
                   do i_b_long = 1, b_num_long
!                      write(nout, *) i_access, i_y, i_b_short, i_b_long
                      i_b_long_global = i_b_long
                      b_long_initial = b_grid_long(i_b_long)
                      i_default_global=2 !Country defaults
                      
                      call optimize_excl(b_next_short, v_valor)
                      b1_next(1) = b_next_short
                      b1_next(2) = b_long_initial*exp(r)
                      disut_default = d0 + d1*y_initial
                      v1_matrix_new(i_b_short, i_b_long, i_y_global, i_access) = v_valor - disut_default
                      b1_next_short_matrix_new(i_b_short, i_b_long, i_y_global, i_access) = b1_next(1)
                      c1_matrix(i_b_short, i_b_long, i_y_global, i_access) = c_global
                      output1_matrix(i_b_short, i_b_long, i_y_global, i_access) = output_global
                       
indicator_tirar = 0

                      i_default_global=1   ! Country does not default and is not excluded for the next period
                      i_excl_global=1
!                      if (i_access ==1) then

                      call optimize(b0_next, v_valor)
                      indicator_tirar = 0
                      v0_matrix_new(i_b_short, i_b_long, i_y, i_access) = v_valor

                       b0_next_short_matrix_new(i_b_short, i_b_long, i_y_global, i_access) = b0_next(1)
                       b0_next_long_matrix_new(i_b_short, i_b_long, i_y_global, i_access)  = b0_next(2)
                       v0_matrix_new(i_b_short, i_b_long, i_y, i_access) = v_valor
                       c0_matrix(i_b_short, i_b_long, i_y_global, i_access) = c_global
                       output0_matrix(i_b_short, i_b_long, i_y_global, i_access) = output_global

                q_nodef_matrix_new(i_b_short, i_b_long, i_y_global, i_access) = q_fun(b0_next(1), b0_next(2), y_initial, i_access)
                      indicator_tirar =0
                 if (v1_matrix_new(i_b_short, i_b_long, i_y, i_access) <= v0_matrix_new(i_b_short, i_b_long, i_y, i_access)) then
                    default_decision(i_b_short, i_b_long, i_y, i_access) = 1
                    b_next_matrix(i_b_short, i_b_long,  i_y, i_access, 1) = b0_next(1)  !SAVINGS IF NO DEFAULT
                    b_next_matrix(i_b_short, i_b_long,  i_y, i_access, 2) = b0_next(2) !SAVINGS IF NO DEFAULT
             q_paid_matrix_new(i_b_short, i_b_long, i_y_global, i_access) =  q_fun(b0_next(1), b0_next(2), y_initial, i_access)
                else
                    default_decision(i_b_short, i_b_long, i_y, i_access) = 2
                    b_next_matrix(i_b_short, i_b_long,  i_y, i_access, 1) = b1_next(1)  !SAVINGS IF NO DEFAULT
                    b_next_matrix(i_b_short, i_b_long,  i_y, i_access, 2) = b1_next(2)  !SAVINGS IF NO DEFAULT
             q_paid_matrix_new(i_b_short, i_b_long, i_y_global, i_access) =  q_fun(b1_next(1), b1_next(2),  y_initial, i_access)
                 end if



                    v_matrix_new(i_b_short, i_b_long, i_y, i_access) = MAX(v1_matrix_new(i_b_short, i_b_long, i_y, i_access), &
                                                                   v0_matrix_new(i_b_short, i_b_long, i_y, i_access))




  deviation = MAX(deviation, ABS(v_matrix_new(i_b_short, i_b_long, i_y, i_access) - v_matrix(i_b_short, i_b_long, i_y, i_access)))
                      dev_matrix(i_b_short, i_b_long, i_y, i_access) = &
                      ABS(v_matrix_new(i_b_short, i_b_long, i_y, i_access) - v_matrix(i_b_short, i_b_long, i_y, i_access))
!                     WRITE(nout, '(F12.8)') dev_matrix(i_b_short, i_b_long, i_y, i_access)
                      indicator_tirar=0
                      dev_matrix_q(i_b_short, i_b_long, i_y, i_access) = min(100d+0, &
ABS(q_paid_matrix_new(i_b_short, i_b_long, i_y, i_access) - q_paid_matrix(i_b_short, i_b_long, i_y, i_access)))
dev_q = MAX(dev_q, ABS(q_paid_matrix_new(i_b_short, i_b_long, i_y, i_access) - q_paid_matrix(i_b_short, i_b_long, i_y, i_access)))
dev_debt = MAX(dev_debt, ABS(b0_next_long_matrix_new(i_b_short, i_b_long, i_y, i_access) - &
                                 b0_next_long_matrix(i_b_short, i_b_long, i_y, i_access)))
dev_res = MAX(dev_res, ABS(b0_next_short_matrix_new(i_b_short, i_b_long, i_y, i_access) - &
                               b0_next_short_matrix(i_b_short, i_b_long, i_y, i_access)))

if(v0_matrix_new(i_b_short, i_b_long, i_y, i_access) < -20000d+0 .or.  v1_matrix_new(i_b_short, i_b_long, i_y, i_access) < -20000d+0) then
pause
end if


                   end do
               end do
         end do
      end do

      
      
      
!close(117)
 WRITE(nout, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F16.8, X, F16.8)') deviation, dev_q 



open(100, FILE ='iteration.txt',POSITION ='append')
write(100, '(F12.8, X, F12.8, X, F12.8, X, F12.8, X, F12.8)') deviation, dev_q, dev_res, dev_debt, counter
close(100)
!3) SAVE RESULTS OF THE CURRENT ITERATION

      open (10, FILE='v.txt',STATUS='replace')
      open (11, FILE='default.txt',STATUS='replace')
      open (12, FILE='q.txt',STATUS='replace')
      open (13, FILE='b_next.txt',STATUS='replace')
      open (14, FILE='dev.txt',STATUS='replace')
      open (16, FILE='q_paid.txt',STATUS='replace')
      open (18, FILE='b_guess.txt',STATUS='replace')
      open (19, FILE='c.txt',STATUS='replace')
      open (23, FILE='counter.txt',STATUS='replace')

      write(23, '(F16.8)') counter
      
      do i_b_short = 1,b_num_short
         do i_b_long = 1,b_num_long
             do i_y = 1,y_num
                do i_access = 1,access_num

                     i_y_global = i_y
                     y_initial = y_grid(i_y)
                     b_short_initial = b_grid_short(i_b_short)
                     i_b_short_global = i_b_short

                     b_long_initial = b_grid_long(i_b_long)
                     i_b_long_global = i_b_long



                     indicator_tirar=0
                     WRITE(10, '(F20.10, X, F20.10, X, F20.10)') v_matrix_new(i_b_short, i_b_long, i_y, i_access), &
                     v0_matrix_new(i_b_short, i_b_long, i_y, i_access), v1_matrix_new(i_b_short, i_b_long, i_y, i_access)
                     WRITE(11, '(F6.2)') default_grid(default_decision(i_b_short, i_b_long, i_y, i_access))

                     d = default_decision(i_b_short, i_b_long, i_y, i_access)

                     g = y_grid(i_y)

                     q_value  = q_fun(b_short_initial, b_long_initial, y_initial, i_access)
                     WRITE(12, '(F15.8, X, F15.8)') q_short_global, q_value
!                   if (i_y == 13 .AND. i_b_short, i_b_long < 40) then
!                      !indicator_tirar = 1
!                   end if

                     WRITE(13, '(F15.11, X, F15.11)') b_next_matrix(i_b_short, i_b_long, i_y, i_access,1), &
                                                      b_next_matrix(i_b_short, i_b_long, i_y, i_access,2)
                     WRITE(14, '(F15.11, X, F15.11)') dev_matrix(i_b_short, i_b_long, i_y, i_access), &
                                                    dev_matrix_q(i_b_short, i_b_long, i_y, i_access)
                 q_value  = q_fun(b_next_matrix(i_b_short, i_b_long, i_y, i_access,1), &
                 b_next_matrix(i_b_short, i_b_long, i_y, i_access,2), y_initial, i_access)
                     WRITE(16, '(F15.8, X, F15.8, X, F15.8)') q_paid_matrix(i_b_short, i_b_long, i_y, i_access),&
                                        q_nodef_matrix_new(i_b_short, i_b_long, i_y, i_access)
                     
                     indicator_tirar=0
 WRITE(18, '(F15.8, X, F15.8, X, F15.8)') b0_next_short_matrix_new(i_b_short, i_b_long, i_y, i_access), &
                                           b0_next_long_matrix_new(i_b_short, i_b_long, i_y, i_access), &
                                          b1_next_short_matrix_new(i_b_short, i_b_long, i_y, i_access)
 
                    WRITE(19, '(F15.8, X, F15.8)') c0_matrix(i_b_short, i_b_long, i_y, i_access), c1_matrix(i_b_short, i_b_long, i_y, i_access)
 
                end do
            end do
         end do
      end do
 CLOSE(10)
 CLOSE(11)
 CLOSE(12)
 CLOSE(13)
 CLOSE(14)
 CLOSE(16)
 close(18)
 close(19)
 close(23)

 

write(nout ,*) 'Iteration saved '
   v0_matrix  = v0_matrix_new
   v1_matrix  = v1_matrix_new
   v_matrix   = v_matrix_new
q_paid_matrix = q_paid_matrix_new
q_nodef_matrix = q_nodef_matrix_new
b0_next_short_matrix = b0_next_short_matrix_new 
b0_next_long_matrix  = b0_next_long_matrix_new
b1_next_short_matrix = b1_next_short_matrix_new

   if (deviation < criteria ) then !.or. counter_final > 15) then
      convergence =1
   end if
end do
end subroutine


program main
USE param
DOUBLE PRECISION :: y_valor, b_valor, f_valor, q_fun, start_time, end_time, indicator_external, def, &
                    b_short, b_long, tax, tfp, labor, output, g, consum, dur_excl, tax0_fun, tax1_fun
INTEGER  i_b_short, i_b_long, i_y, i_access, i
EXTERNAL q_fun


open (1987, FILE='sigma.txt')
read(1987, '(F10.4)') sigma
close(1987)

open (1987, FILE='beta.txt')
read(1987, '(F10.4)') beta
close(1987)

open (1987, FILE='d0.txt')
read(1987, '(F10.4)') d0
close(1987)

open (1987, FILE='d1.txt')
read(1987, '(F10.4)') d1
close(1987)

open (1987, FILE='dur_excl.txt')
read(1987, '(F10.4)') dur_excl
close(1987)

open (1987, FILE='delta.txt')
read(1987, '(F10.4)') delta
close(1987)

open (1987, FILE='premium.txt')
read(1987, '(F10.4)') premia_grid(2)
close(1987)

open (1987, FILE='output_cost.txt')
read(1987, '(F10.4)') output_cost
close(1987)


!d0 = d0-d1
b_recov = 0d+0 !zero recovery
prob_excl_end = 1d+0/ dur_excl

trans_matrix(2,2) = 0.2d+0 !PROBABILITY OF REMAINING IN HIGH PREMIA STATE
trans_matrix(2,1) = 1d+0 - trans_matrix(2,2)  !PROBABILITY OF SWITCHING FORM HIGH PREMIA to LOW PREMIA
trans_matrix(1,1) = 0.85d+0 !    !PROBABILITY OF REMAINING IN LOW PREMIA STATE
trans_matrix(1,2) = 1d+0 - trans_matrix(1,1)    !PROBABILITY OF SWITCHING FORM LOW PREMIA to HIGH PREMIA

premia_grid(1) = 0d+0



write(6, '(F17.4)') sigma
write(6, '(F17.4)') beta
write(6, '(F17.4)') d0
write(6, '(F17.4)') d1
write(6, '(F17.4)') prob_excl_end
write(6, '(F17.4)') delta
write(6, '(F17.4)') premia_grid(2)
write(6, '(F17.4)') output_cost


 
call cpu_time(start_time)
call compute_grid
call quadrature
call compute_knots



indicator_external = 0 !FROM EXTERNAL FILE

    
  if (indicator_external < 0.5) then
        counter = 0
        
         do i_access = 1,access_num
            i_access_initial = i_access
            do i_y = 1,y_num
                y_initial = y_grid(i_y)

                do i_b_short = 1,b_num_short
                   b_short_initial = b_grid_short(i_b_short)
                   b_short = b_inf_short
                   b_long = b_inf
                   i_default_global = 2
                   
                   output = exp(y_initial)*(1-output_cost)  
                   consum = output + b_short_initial - absortion
                   disut_default = d0 + d1*y_initial
                   v1_matrix(i_b_short, :, i_y, i_access) = (consum)**(1-sigma)/(1-sigma) - disut_default 
 !           WRITE(nout, '(F14.8, X, F14.8, X, F14.8, X, F14.8, X, F14.8, X, F14.8)') g_initial-b_short_initial, tax, tfp, output, consum, v1_matrix(i_b_short, 1, i_y, i_access)

                   do i_b_long = 1,b_num_long
                      i_default_global = 1
                      b_long_initial = b_grid_long(i_b_long)
                      output = EXP(y_initial) 
                      consum = output + b_short_initial - coupon*b_long_initial - absortion
                      v0_matrix(i_b_short, i_b_long, i_y, i_access) = (consum)**(1-sigma)/(1-sigma)
!WRITE(nout, '(F14.8, X, F14.8, X, F14.8, X, F14.8, X, F14.8, X, F14.8)') output, consum, v0_matrix(i_b_short, i_b_long, i_y, i_access), (consum)**(1-sigma)/(1-sigma)
                      
                      v_matrix(i_b_short, i_b_long, i_y, i_access) = MAX(v1_matrix(i_b_short, i_b_long,i_y, i_access), &
                                                                         v0_matrix(i_b_short, i_b_long, i_y, i_access))
                     
                   end do
                end do
            end do
         end do
   else !READ DATA FROM EXTERNAL FILES

       open (10, FILE='v.txt')
       open (16, FILE='q_paid.txt')
       open (18, FILE='b_guess.txt')
       open (19, FILE='counter.txt')

       read (19, '(F16.8)') counter

     do i_b_short = 1,b_num_short
         do i_b_long = 1,b_num_long
             do i_y = 1,y_num
                do i_access = 1,access_num
                READ(10, '(F20.10, X, F20.10, X, F20.10)') v_matrix(i_b_short, i_b_long, i_y, i_access), &
                    v0_matrix(i_b_short, i_b_long, i_y, i_access), v1_matrix(i_b_short, i_b_long, i_y, i_access)

                if (v0_matrix(i_b_short, i_b_long, i_y, i_access)< v1_matrix(i_b_short, i_b_long, i_y, i_access)) then
                   default_decision(i_b_short, i_b_long, i_y, i_access) =2
                else
                   default_decision(i_b_short, i_b_long, i_y, i_access) =1
                end if
                READ(16, '(F15.8, X, F15.8, X, F15.8)')  q_paid_matrix(i_b_short, i_b_long, i_y, i_access), &
                                           q_nodef_matrix(i_b_short, i_b_long, i_y, i_access)
READ(18, '(F15.8, X, F15.8, X, F15.8)') b0_next_short_matrix(i_b_short, i_b_long, i_y, i_access), &
                                         b0_next_long_matrix(i_b_short, i_b_long, i_y, i_access), &
                                        b1_next_short_matrix(i_b_short, i_b_long, i_y, i_access)

 
                end do
              end do
         end do
     end do
     CLOSE(10)
     CLOSE(16)
     close(18)
     close(19)

end if



call iterate

call cpu_time(end_time)
WRITE(nout, '(A7, X, A7, X, A7)') 'Hours ', 'Minutes', 'Seconds'
WRITE(nout, '(I7, X, I7, X, I7)') INT((end_time - start_time) / 3600d+0), &
             INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0),&
INT(end_time-start_time - INT((end_time - start_time) / 3600d+0)*3600d+0 - &
INT((end_time-start_time)/60d+0 - INT((end_time - start_time) / 3600d+0)*60d+0)*60d+0)
end program




SUBROUTINE POWELL(FUNC,P,XI,N,NP,FTOL,ITER,FRET)
!-----------------------------------------------------------
! Minimization of a function  FUNC of N variables  (FUNC is
! not an argument, it is a fixed function name). Input con-
! sists of an initial starting point P  that is a vector of
! length N; an initial matrix XI  whose  logical dimensions
! are N by N, physical dimensions NP by NP, and whose columns26
! contain the initial set of directions (usually the N unit
! vectors); and FTOL, the fractional tolerance in the func-
! tion value such that failure to decrease by more than this
! amount on one iteration signals doneness. On output, P is
! set to the best point found, XI is the then-current direc-
! tion set,  FRET is the returned function value at P,  and
! ITER is the number of iterations taken. The routine LINMIN
! is used.
!------------------------------------------------------------
  !IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
  INTEGER :: N, NP, NMAX, ITMAX, I, IBIG, J, ITER
  PARAMETER(NMAX=2,ITMAX=100)
  DOUBLE PRECISION :: P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)
  DOUBLE PRECISION :: FRET, FP, DEL, FPTT, FTOL, T
  !EXTERNAL objective_function_nobind
  
  interface
  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
    DOUBLE PRECISION, intent(in) :: ARG(2)
  end function
  end interface
  
  FRET=  FUNC(P) !objective_function_nobind(P) !
  !write(nout, *) FRET, objective_function_nobind(P)
  !pause
  
  DO J=1,N
    PT(J)=P(J)       !Save initial pont
  END DO
  ITER=0
1 ITER=ITER+1
  FP=FRET
  IBIG=0
  DEL=0.0D+0           !Will be the biggest function decrease.
  DO I=1,N           !In each iteration, loop over all directions in the set.
    DO J=1,N         !Copy the direction
      XIT(J)=XI(J,I)
    END DO
    FPTT=FRET
    CALL LINMIN(FUNC,P,XIT,N,FRET)  !Minimize along it.
    IF (DABS(FPTT-FRET).GT.DEL) THEN
      DEL=DABS(FPTT-FRET)
      IBIG=I
    END IF
  END DO
  IF (2.D0*DABS(FP-FRET).LE.FTOL*(DABS(FP)+DABS(FRET))) RETURN !Termination criterion
  IF (ITER.EQ.ITMAX) Then
    Pause ' Powell exceeding maximum iterations.'
    return
  END IF		 
  DO J=1,N
    PTT(J)=2.D0*P(J)-PT(J)  !Construct the extrapolated point and the average
    XIT(J)=P(J)-PT(J)       !direction moved. Save the old starting point.
    PT(J)=P(J)
  END DO
  FPTT= FUNC(PTT) ! objective_function_nobind(PTT) !          !Function value at extrapolated point.
!  write(nout, '(F12.8, X, F12.8, X, F15.6, X, F15.6)') PTT(1), PTT(2), FPTT, objective_function_nobind(PTT)
  IF (FPTT.GE.FP) GO TO 1   !One reason not to use new direction. 
  T=2.D0*(FP-2.D0*FRET+FPTT)*(FP-FRET-DEL)**2D+0-DEL*(FP-FPTT)**2D+0
  IF (T.GE.0.D0) GO TO 1    !Other reason not to use new direction.
  CALL LINMIN(FUNC,P,XIT,N,FRET) !Move to the minimum of the new direction.
  DO J=1,N                  !and save the new direction
    XI(J,IBIG)=XIT(J)
  END DO
  GO TO 1
END

SUBROUTINE LINMIN(FUNC,P,XI,N,FRET)
!----------------------------------------------------------
! Given an N dimensional point P and a N dimensional direc-
! tion XI, moves and resets P to where the function FUNC(P)
! takes on a minimum along the direction XI from P, and 
! replaces XI by the actual vector displacement that P was
! moved. Also returns as FRET the value of FUNC at the
! returned location P. This is actually all accomplished by
! calling the routines MNBRAK and BRENT.
!----------------------------------------------------------
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DOUBLE PRECISION :: TOL
  INTEGER :: N, NMAX, NCOM, J
  PARAMETER(NMAX=50,TOL=1.D-4)
  DOUBLE PRECISION :: P(N),XI(N) !, PCOM(N), XICOM(N)
  DOUBLE PRECISION :: AX, XX, BX, FA, FX,FB, FRET, BRENT
  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM

  interface
  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
    DOUBLE PRECISION, intent(in) :: ARG(2)
  end function
  end interface
  NCOM=N
  DO J=1,N
    PCOM(J)=P(J)
    XICOM(J)=XI(J)
  END DO
  AX=0.D0
  XX=1.D0
  BX=2.D0
  CALL MNBRAK(FUNC,AX,XX,BX,FA,FX,FB)
  FRET=BRENT(FUNC,AX,XX,BX,TOL,XMIN)
  DO J=1,N
    XI(J)=XMIN*XI(J)
    P(J)=P(J)+XI(J)
  END DO
  RETURN
END

DOUBLE PRECISION FUNCTION F1DIM(FUNC,X)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER :: NMAX, J, NCOM
  PARAMETER(NMAX=50)
  COMMON /F1COM/ PCOM(NMAX),XICOM(NMAX),NCOM
  DOUBLE PRECISION :: XT(NMAX), X
  !EXTERNAL objective_function_nobind
 ! EXTERNAL FUNC

  interface
  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
    DOUBLE PRECISION, intent(in) :: ARG(2)
  end function
  end interface
  
  DO J=1, NCOM
    XT(J)=PCOM(J)+X*XICOM(J)
  END DO
  F1DIM = FUNC(XT) !objective_function_nobind(XT) !
  RETURN
END

SUBROUTINE MNBRAK(FUNC,AX,BX,CX,FA,FB,FC)
!----------------------------------------------------------------------
!Given a Function F1DIM(X), and given distinct initial points AX and
!BX, this routine searches in the downhill direction (defined by the
!F1DIMtion as evaluated at the initial points) and returns new points
!AX, BX, CX which bracket a minimum of the Function. Also returned
!are the Function values at the three points, FA, FB and FC.
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DOUBLE PRECISION :: GOLD, GLIMIT, TINY, FA, FB,FC, AX,BX,CX, DUM,R,Q,U,ULIM,FU, F1DIM
PARAMETER(GOLD=1.618034,GLIMIT=100.,TINY=1.D-20)
  interface
  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
    DOUBLE PRECISION, intent(in) :: ARG(2)
  end function
  end interface

!The first parameter is the default ratio by which successive intervals
!are magnified; the second is the maximum magnification allowed for
!a parabolic-fit step.
!----------------------------------------------------------------------
FA=F1DIM(FUNC,AX)
FB=F1DIM(FUNC,BX)
IF(FB.GT.FA) THEN
  DUM=AX
  AX=BX
  BX=DUM
  DUM=FB
  FB=FA
  FA=DUM
ENDIF
CX=BX+GOLD*(BX-AX)
FC=F1DIM(FUNC,CX)
1 IF(FB.GE.FC) THEN
  R=(BX-AX)*(FB-FC)
  Q=(BX-CX)*(FB-FA)
  U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
  ULIM=BX+GLIMIT*(CX-BX)
  IF((BX-U)*(U-CX).GT.0) THEN
    FU=F1DIM(FUNC,U)
    IF(FU.LT.FC) THEN
      AX=BX
      FA=FB
      BX=U
      FB=FU
      GOTO 1
    ELSE IF(FU.GT.FB) THEN
      CX=U
      FC=FU
      GOTO 1
    ENDIF
    U=CX+GOLD*(CX-BX)
    FU=F1DIM(FUNC,U)
  ELSE IF((CX-U)*(U-ULIM).GT.0) THEN
    FU=F1DIM(FUNC,U)
    IF(FU.LT.FC) THEN
      BX=CX
      CX=U
      U=CX+GOLD*(CX-BX)
      FB=FC
      FC=FU
      FU=F1DIM(FUNC,U)
    ENDIF
  ELSE IF((U-ULIM)*(ULIM-CX).GE.0) THEN
    U=ULIM
    FU=F1DIM(FUNC,U)
  ELSE
    U=CX+GOLD*(CX-BX)
    FU=F1DIM(FUNC,U)
  ENDIF
  AX=BX
  BX=CX
  CX=U
  FA=FB
  FB=FC
  FC=FU
  GOTO 1
ENDIF
RETURN
END

DOUBLE PRECISION FUNCTION BRENT(FUNC,AX,BX,CX,TOL,XMIN)
!-------------------------------------------------------------------
!Given a function F1DIM, and a bracketing triplet of abscissas
!AX,BX,CX (such that BX is between AX and CX, and F(BX) is less 
!than both F(AX) and F(CX)), this routine isolates the minimum 
!to a fractional precision of about TOL using Brent's method.
!The abscissa of the minimum is returned in XMIN, and the minimum
!function value is returned as BRENT, the returned function value.
!-------------------------------------------------------------------
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
INTEGER :: ITMAX,ITER
DOUBLE PRECISION :: CGOLD, ZEPS
PARAMETER(ITMAX=100,CGOLD=.3819660,ZEPS=1.D-10)
!Maximum allowed number of iterations; golden ration; and a small
!number which protects against trying to achieve fractional accuracy
!for a minimum that happens to be exactly zero.
DOUBLE PRECISION :: AX,BX, CX,TOL,XMIN,A,B,V,W,X,E,FX,FV,FU,FW,XM,TOL1,TOL2,R,P,Q,ETEMP,D,U
interface
  DOUBLE PRECISION PURE FUNCTION FUNC(ARG)
    DOUBLE PRECISION, intent(in) :: ARG(2)
  end function
  end interface


A=MIN(AX,CX)
B=MAX(AX,CX)
V=BX
W=V
X=V
E=0.
FX=F1DIM(FUNC,X)
FV=FX
FW=FX
DO 11 ITER=1,ITMAX	                                !main loop
  XM=0.5*(A+B)
  TOL1=TOL*ABS(X)+ZEPS
  TOL2=2.*TOL1
  IF (ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3  !Test for done here
  IF (ABS(E).GT.TOL1) THEN     !Construct a trial parabolic fit
    R=(X-W)*(FX-FV)
    Q=(X-V)*(FX-FW)
    P=(X-V)*Q-(X-W)*R
    Q=.2*(Q-R)
    IF (Q.GT.0)  P=-P
    Q=ABS(Q)
    ETEMP=E
    E=D
    IF (ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.  &
	P.GE.Q*(B-X))  GOTO 1
!   The above conditions determine the acceptability of the 
!   parabolic fit. Here it is o.k.:
    D=P/Q
    U=X+D
    IF (U-A.LT.TOL2.OR.B-U.LT.TOL2)  D=SIGN(TOL1,XM-X)
    GOTO 2
  ENDIF
1 IF (X.GE.XM) THEN
    E=A-X
  ELSE
    E=B-X
  ENDIF
  D=CGOLD*E
2 IF (ABS(D).GE.TOL1) THEN
    U=X+D
  ELSE
    U=X+SIGN(TOL1,D)
  ENDIF
  FU=F1DIM(FUNC,U)  !This the one function evaluation per iteration
  IF (FU.LE.FX) THEN
    IF (U.GE.X) THEN
      A=X
    ELSE
      B=X
    ENDIF
    V=W
    FV=FW
    W=X
    FW=FX
    X=U
    FX=FU
  ELSE
    IF (U.LT.X) THEN
      A=U
    ELSE
      B=U
    ENDIF
    IF (FU.LE.FW.OR.W.EQ.X) THEN
      V=W
      FV=FW
      W=U
      FW=FU
    ELSE IF (FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
      V=U
      FV=FU
    ENDIF
  ENDIF
11 CONTINUE
!Pause ' Brent exceed maximum iterations.'
3 XMIN=X   !exit section
  BRENT=FX
  RETURN
  END
!end of file tpowell.f90

subroutine BISEC(func, x1, x2, ERRABS, ERRREL)
implicit none
DOUBLE PRECISION, INTENT(IN) ::  ERRABS, ERRREL
DOUBLE PRECISION, INTENT(IN OUT) :: x2, x1
DOUBLE PRECISION :: diff
INTEGER :: nout

interface
        function func(x)
        implicit none
        DOUBLE PRECISION, INTENT(IN) :: x
        DOUBLE PRECISION :: func
        end function func
END interface
DOUBLE PRECISION :: fl, fh, f_mean, xl, xh, xmean, rel_error

xl=x1
xh=x2
fl=func(xl)
fh=func(xh)
diff = fh-fl
rel_error=1

do WHILE( MIN((ABS(fh)-ERRABS), rel_error-ERRREL) > 0)

xmean =  0.5d+0 * (xl+xh)
!xl - fl / ((fh - fl)/(xh-xl) )  !Use linear approx to find the zero !(xl+xh)/2
f_mean = func(xmean)

if (fl*fh>0) then
   !WRITE(nout,*) 'Wrong range'       !Need change of signs to compute the root
   !call BREAK
else
   if (diff < 0) then
      if (f_mean > 0) then
         xl=xmean
         fl=f_mean
      else
         xh=xmean
         fh=f_mean
      end if
   else
      if (f_mean > 0) then
         xh=xmean
         fh=f_mean
      else
         xl=xmean
         fl=f_mean
      end if
    end if

end if
rel_error=ABS(xh-xl)

!WRITE(6,*) rel_error, fh
end do
x1 = xl
x2 = xh


end subroutine    
    


    