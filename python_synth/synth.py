import numpy as np
import scipy

data = np.load("data.npz")
alpha, beta, envelope = data['alpha'], data['beta'], data['envelope']

fs = 44100
gamma = 40000.0

ovfs=20
prct_noise=0
t, tmax, dt = 0, int(alpha.size)*ovfs-1, 1./(ovfs*fs)
# pback and pin vectors initialization
pi, pb, out = np.zeros(tmax), np.zeros(tmax), np.zeros(int(alpha.size))
# initial vector ODEs (v0), it is not too relevant
v = 1e-4*np.array([1e2, 1e1, 1, 1, 1, 1]);
# ------------- BIRD PARAMETERS -----------
c, L, r, Ch, MG, MB, RB, Rh = 343.0, 0.025, 0.65, 1.43e-10, 20.0, 10000.0, 5000000.0, 24000.0

# - Trachea:
#           r: reflection coeficient    [adimensionelss]
#           L: trachea length           [m]
#           c: speed of sound in media  [m/s]
# - Beak, Glottis and OEC:
#           CH: OEC Compliance          [m^3/Pa]
#           MB: Beak Inertance          [Pa s^2/m^3 = kg/m^4]
#           MG: Glottis Inertance       [Pa s^2/m^3 = kg/m^4]
#           RB: Beak Resistance         [Pa s/m^3 = kg/m^4 s]
#           Rh: OEC Resistence          [Pa s/m^3 = kg/m^4 s]
# ------------------------------ ODEs -----------------------------
def ODEs(v):
    dv, [x, y, pout, i1, i2, i3] = np.zeros(6), v  # (x, y, pout, i1, i2, i3)'
    # ----------------- direct implementation of the EDOs -----------
    dv[0] = y
    dv[1] = (-alpha[t//ovfs]-beta[t//ovfs]*x-x**3+x**2)*gamma**2 - (x**2*y+x*y)*gamma
    # ------------------------- trachea ------------------------
    pbold = pb[t]                                 # pressure back before
    # Pin(t) = Ay(t)+pback(t-L/C) = envelope_Signal*v[1]+pb[t-L/C/dt]
    pi[t] = (.5*envelope[t//ovfs])*dv[1] + pb[t-int(L/c/dt)] 
    pb[t] = -r*pi[t-int(L/c/dt)]                          # pressure back after: -rPin(t-L/C) 
    pout  = (1-r)*pi[t-int(L/c/dt)]                       # pout
    # ---------------------------------------------------------------
    dv[2] = (pb[t]-pbold)/dt                      # dpout
    dv[3] = i2
    dv[4] = -(1/Ch/MG)*i1 - Rh*(1/MB+1/MG)*i2 +(1/MG/Ch+Rh*RB/MG/MB)*i3 \
            +(1/MG)*dv[2] + (Rh*RB/MG/MB)*pout
    dv[5] = -(MG/MB)*i2 - (Rh/MB)*i3 + (1/MB)*pout
    return dv

def rk4(f, v, dt):
    """
    Implent of Runge-Kuta 4th order
    INPUT:
        f  = differential equations functions y'=f(y)
        v  = vector y of the differential equations [x,y,i1,i2,i3]
        dt = rk4 time step
    OUTPUT:
        rk4 approximation 
    """
    k1 = f(v)    
    k2 = f(v + dt/2.0*k1)
    k3 = f(v + dt/2.0*k2)
    k4 = f(v + dt*k3)
    return v + dt*(2.0*(k2+k3)+k1+k4)/6.0

# ----------------------- Solving EDOs ----------------------
while t < tmax: # and v[1] > -5e6:  # labia velocity not too fast
    v = rk4(ODEs, v, dt);  # self.Vs.append(v)  # RK4 - step
    out[t//ovfs] = RB*v[-1]               # output signal (synthetic) 
    t += 1;
    
    # # if the bird OEC change of size in time
    # BirdData = pd.read_csv(self.paths.auxdata/'ZonotrichiaData.csv')
    # c, L, r, Ch, MG, MB, RB, Rh = BirdData['value'] # c, L, r, c, L1, L2, r2, rd 

# ------------------------------------------------------------

print(out)
scipy.io.wavfile.write("out.wav", fs, out / np.max(np.abs(out)))
