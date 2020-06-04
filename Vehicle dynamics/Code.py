# for loaded vehicle
# Given parameters

w1 = 25600                                      # weight of vehicle on axle 1 (in N)
w2 = 33600                                      # weight of vehicle on axle 2 (in N)
k_s1 = 79200                                    # front axle sprung stiffness(N/m)
k_s2 = 137200                                   # rear axle sprung stiffness(N/m)
k_t = 1618000                                   # Vertical stiffness of tire for front and rear axle (N/m)
h1 = 2                                          # center of mass height of front axle
h2 = 2.6                                        # center of mass height of rear axle
t =1                                            # in m
a_y_def = 0
h21= 0.65                                       # roll center height of front axle
h22 = 0.65                                      # roll center height of rear axle

def theta_s(w,t,k_s,k_t,a_y):
    return ((w*t)*((1/k_s)+(1/k_t))) - ((w*a_y*h21)/k_s)                # roll angle of sprung mass formula

def moment(w1,w2,t1,t2,h1,h2,theta_s1,a_y):
  return ((w1*t1+w2*t2) - ((w1*h1+w2*h2)*theta_s1))/ (w1*h1+w2*h2)      # moment equation

flag = 1

while(flag ==1):

  theta_s1 = theta_s(w1,t,k_s1,k_t,a_y_def)                             # roll angle of sprung mass for front axle
  theta_s2 = theta_s(w2,t,k_s2,k_t,a_y_def)                             # roll angle of sprung mass for rear axle

  if theta_s1> theta_s2:
    a_y_new = round(moment(w1,w2,t,t,h1,h2,theta_s1,a_y_def),3)
    if a_y_new == a_y_def:
      flag = 0
    else:
      a_y_def = round(a_y_def+0.001,3)
      if(a_y_def>1):
        flag=0

print("theta_s1={0}".format(theta_s1))
print("theta_s2={0}".format(theta_s2))
print("a_y={0}".format(a_y_new))

