#include <stdio.h>
#include <math.h>

static const double GRAVITY = 32.2; // ft/s^2
static const double PHI = 1.486;

static double area(double width, double y){
    return width * y;
}
static double hyd_radius(double width, double y){
    return width*y/(width + 2.0*y);
}

static double compute_step(double q_last, double y1, double y2, double h1, double h2,
                           double q_old, double a_old, double dt, double length,
                           double width, double n, double omega){
    double a1 = area(width, y1);
    double a2 = area(width, y2);
    double a_mid = 0.5*(a1+a2);
    double r1 = hyd_radius(width, y1);
    double r_mid = hyd_radius(width, 0.5*(y1+y2));

    double rho = 1.0; // no damping selection
    double a_wtd = a1 + (a_mid - a1)*rho;
    double r_wtd = r1 + (r_mid - r1)*rho;

    double v = q_last / a_mid;
    double roughFactor = GRAVITY * (n/PHI)*(n/PHI);

    double dq1 = dt * roughFactor / pow(r_wtd, 1.33333) * fabs(v);
    double dq2 = dt * GRAVITY * a_wtd * (h2 - h1) / length;
    double denom = 1.0 + dq1;
    double q = (q_old - dq2) / denom;

    // under-relaxation
    q = (1.0 - omega)*q_last + omega*q;
    return q;
}

int main(){
    double width=3.0, length=200.0, n=0.013, dt=1.0, omega=0.5;
    double y1=5.0, y2=4.9, h1=5.0, h2=4.9; // depths & heads above invert
    double q_old=0.0, a_old=area(width, (y1+y2)/2.0);

    double q_last = q_old;
    for(int iter=0; iter<20; iter++){
        double q_new = compute_step(q_last,y1,y2,h1,h2,q_old,a_old,dt,length,width,n,omega);
        if (fabs(q_new - q_last) < 1e-6) break;
        q_last = q_new;
    }
    printf("%f\n", q_last);
    return 0;
}
