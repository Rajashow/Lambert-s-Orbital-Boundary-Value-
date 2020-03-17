def solve_lambert(r1,r2,t,mu):
  if not (t>0 and mu>0):
  #     throw error
      pass
  c = diff(r1,r2)
  c = norm_abs(c)
  r1_ = norm_abs(r1)
  r2_ = norm_abs(r2)
  s = sum_(r1_,r2_,c)/2
  i_r_1,i_r_2 = r1/r1_, r2/r2_
  i_h = mult_(i_r_1,i_r_2)
  lambda = sqrt(1- (c/s))
  if diff(mult(r1[0],r1[1]),mult(r2[0],r2[1]) ) <0:
    lambda = -lambda
    i_t_1 = mult(i_r_1,i_h)
    i_t_2 = mult(i_r_2,i_r_2)
  else:
    i_t_1 = mult(i_r_1,i_h)
    i_t_2 = mult(i_r_2,i_r_2)
  T = sqrt(2*mu/ cube(s))*t
  xs,ys = findxy(lambda,T)
  gamma = sqrt(mu*s/2)
  rho = diff(r1,r2)/c
  sigma = sqrt(1-rho**2)
  for x,y in zip(xs,ys):
    num = gamma* (((lambda*y) -x)- rho*(lambda*y + x))
    v_r_1 = num/r1_
    v_r_2 = -num/r2_
    vt_num = (y+ lambda*x)*gamma*sigma
    v_t1 = vt_num/r1
    v_t2 = vt_num/r2
    v1 = v_r_1*i_r_1 + v_t1*i_t_1
    v2 = v_r_2*i_r_2 + v_t2*i_t_2
  pass
