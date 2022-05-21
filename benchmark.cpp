#include "WignerSymbol.hpp"
#include <chrono>
#include <gsl/gsl_specfunc.h>

using namespace jshl;

void time_6j() {
  using timer_colok = std::chrono::high_resolution_clock;
  auto t1 = timer_colok::now();
  WignerSymbols wigner;
  int N = 20;
  wigner.reserve(N, "Jmax", 6);
  double x = 0;
  for (int dj1 = N; dj1 <= 2*N; ++dj1) {
    for (int dj2 = N; dj2 <= 2*N; ++dj2) {
      for (int dj3 = N; dj3 <= 2*N; ++dj3) {
        for (int dj4 = N; dj4 <= 2*N; ++dj4) {
          for (int dj5 = N; dj5 <= 2*N; ++dj5) {
            for (int dj6 = N; dj6 <= 2*N; ++dj6) {
              x += wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
            }
          }
        }
      }
    }
  }
  auto t2 = timer_colok::now();
  double y = 0;
  for (int dj1 = N; dj1 <= 2*N; ++dj1) {
    for (int dj2 = N; dj2 <= 2*N; ++dj2) {
      for (int dj3 = N; dj3 <= 2*N; ++dj3) {
        for (int dj4 = N; dj4 <= 2*N; ++dj4) {
          for (int dj5 = N; dj5 <= 2*N; ++dj5) {
            for (int dj6 = N; dj6 <= 2*N; ++dj6) {
              y += gsl_sf_coupling_6j(dj1, dj2, dj3, dj4, dj5, dj6);
            }
          }
        }
      }
    }
  }
  auto t3 = timer_colok::now();
  std::cout << "time 6j, diff = " << x - y << std::endl;
  std::cout
      << "this code time = "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
      << " ms" << std::endl;
  std::cout
      << "gsl library time = "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
      << " ms" << std::endl;
}

void time_9j() {
  using timer_colok = std::chrono::high_resolution_clock;
  auto t1 = timer_colok::now();
  WignerSymbols wigner;
  int N = 8;
  wigner.reserve(N, "Jmax", 9);
  double x = 0;
  for (int dj1 = N; dj1 <= 2*N; ++dj1) {
    for (int dj2 = N; dj2 <= 2*N; ++dj2) {
      for (int dj3 = N; dj3 <= 2*N; ++dj3) {
        for (int dj4 = N; dj4 <= 2*N; ++dj4) {
          for (int dj5 = N; dj5 <= 2*N; ++dj5) {
            for (int dj6 = N; dj6 <= 2*N; ++dj6) {
              for (int dj7 = N; dj7 <= 2*N; ++dj7) {
                for (int dj8 = N; dj8 <= 2*N; ++dj8) {
                  for (int dj9 = N; dj9 <= 2*N; ++dj9) {
                    x +=
                        wigner.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  auto t2 = timer_colok::now();
  double y = 0;
  for (int dj1 = N; dj1 <= 2*N; ++dj1) {
    for (int dj2 = N; dj2 <= 2*N; ++dj2) {
      for (int dj3 = N; dj3 <= 2*N; ++dj3) {
        for (int dj4 = N; dj4 <= 2*N; ++dj4) {
          for (int dj5 = N; dj5 <= 2*N; ++dj5) {
            for (int dj6 = N; dj6 <= 2*N; ++dj6) {
              for (int dj7 = N; dj7 <= 2*N; ++dj7) {
                for (int dj8 = N; dj8 <= 2*N; ++dj8) {
                  for (int dj9 = N; dj9 <= 2*N; ++dj9) {
                    y += gsl_sf_coupling_9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7,
                                            dj8, dj9);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  auto t3 = timer_colok::now();
  std::cout << "time 9j, diff = " << x - y << std::endl;
  std::cout
      << "this code time = "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count()
      << " ms" << std::endl;
  std::cout
      << "gsl library time = "
      << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()
      << " ms" << std::endl;
}

void test_3j() {
  WignerSymbols wigner;
  int N = 10;
  wigner.reserve(N, "Jmax", 3);
  double diff = 0;
  for (int dj1 = N; dj1 <= 2*N; ++dj1) {
    for (int dj2 = N; dj2 <= 2*N; ++dj2) {
      for (int dj3 = N; dj3 <= 2*N; ++dj3) {
        for (int dm1 = -dj1; dm1 <= dj1; ++dm1) {
          for (int dm2 = -dj2; dm2 <= dj2; ++dm2) {
            for (int dm3 = -dj3; dm3 <= dj3; ++dm3) {
              double x = wigner.f3j(dj1, dj2, dj3, dm1, dm2, dm3);
              double y = gsl_sf_coupling_3j(dj1, dj2, dj3, dm1, dm2, dm3);
              diff += std::abs(x - y);
            }
          }
        }
      }
    }
  }
  std::cout << "test 3j, diff = " << diff << std::endl;
}

void test_6j() {
  WignerSymbols wigner;
  int N = 20;
  wigner.reserve(N, "Jmax", 6);
  double diff = 0;
  for (int dj1 = N; dj1 <= 2*N; ++dj1) {
    for (int dj2 = N; dj2 <= 2*N; ++dj2) {
      for (int dj3 = N; dj3 <= 2*N; ++dj3) {
        for (int dj4 = N; dj4 <= 2*N; ++dj4) {
          for (int dj5 = N; dj5 <= 2*N; ++dj5) {
            for (int dj6 = N; dj6 <= 2*N; ++dj6) {
              double x = wigner.f6j(dj1, dj2, dj3, dj4, dj5, dj6);
              double y = gsl_sf_coupling_6j(dj1, dj2, dj3, dj4, dj5, dj6);
              diff += std::abs(x - y);
            }
          }
        }
      }
    }
  }
  std::cout << "test 6j, diff = " << diff << std::endl;
}

void test_9j() {
  WignerSymbols wigner;
  int N = 6;
  wigner.reserve(N, "Jmax", 9);
  double diff = 0;
  for (int dj1 = 0; dj1 <= N; ++dj1) {
    for (int dj2 = 0; dj2 <= N; ++dj2) {
      for (int dj3 = 0; dj3 <= N; ++dj3) {
        for (int dj4 = 0; dj4 <= N; ++dj4) {
          for (int dj5 = 0; dj5 <= N; ++dj5) {
            for (int dj6 = 0; dj6 <= N; ++dj6) {
              for (int dj7 = 0; dj7 <= N; ++dj7) {
                for (int dj8 = 0; dj8 <= N; ++dj8) {
                  for (int dj9 = 0; dj9 <= N; ++dj9) {
                    double x =
                        wigner.f9j(dj1, dj2, dj3, dj4, dj5, dj6, dj7, dj8, dj9);
                    double y = gsl_sf_coupling_9j(dj1, dj2, dj3, dj4, dj5, dj6,
                                                  dj7, dj8, dj9);
                    diff += std::abs(x - y);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  std::cout << "test 9j, diff = " << diff << std::endl;
}

int main() {
  test_3j();
  test_6j();
  test_9j();
  time_6j();
  time_9j();
  return 0;
}
