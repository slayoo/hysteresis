#pragma once

// identity: xi = rw
template <typename T=double> struct xi_id
{
  static T xi_of_rw1(const T &rw) { return rw; }
  static T xi_of_rw3(const T &rw3) { return pow(rw3, T(1./3)); }
  static T rw1_of_xi(const T &xi) { return xi; }
  static T rw2_of_xi(const T &xi) { return xi * xi; }
  static T rw3_of_xi(const T &xi) { return xi * xi * xi; }
  static T dxidrw(const T &) { return 1; }
};

// xi = ln(rw / nm)
template <typename T=double> struct xi_ln
{
  static T xi_of_rw1(const T &rw) { return log(rw / T(1e-9)); }
  static T xi_of_rw3(const T &rw3) { return log(pow(rw3, T(1./3)) / T(1e-9)); }
  static T rw1_of_xi(const T &xi) { return T(1e-9) * exp(xi); }
  static T rw2_of_xi(const T &xi) { return pow(T(1e-9) * exp(xi), 2); }
  static T rw3_of_xi(const T &xi) { return pow(T(1e-9) * exp(xi), 3); }
  static T dxidrw(const T &rw) { return T(1) / rw; }
};

// xi = pow(rw, 2)
template <typename T=double> struct xi_p2
{
  static T xi_of_rw1(const T &rw) { return rw * rw; }
  static T xi_of_rw3(const T &rw3) { return pow(rw3, T(2./3)); }
  static T rw1_of_xi(const T &xi) { return sqrt(xi); }
  static T rw2_of_xi(const T &xi) { return xi; }
  static T rw3_of_xi(const T &xi) { return pow(xi, T(3./2)); }
  static T dxidrw(const T &rw) { return T(2) * rw; }
};

// xi = pow(rw, 3)
template <typename T=double> struct xi_p3
{
  static T xi_of_rw1(const T &rw) { return rw * rw * rw; }
  static T xi_of_rw3(const T &rw3) { return rw3; }
  static T rw1_of_xi(const T &xi) { return pow(xi, T(1)/T(3)); }
  static T rw2_of_xi(const T &xi) { return pow(xi, 2./3); }
  static T rw3_of_xi(const T &xi) { return xi; }
  static T dxidrw(const T &rw) { return T(3) * rw * rw; }
};

