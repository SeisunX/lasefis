#include <rsf.h>
#include <assert.h>


/**
 * Addressing the array in a plane-1D way. It is CPE-friendly.
 * We need to make sure that the size of the array is (nx2 * ny2 * nz2),
 */
void fd2_addr1d(float DT, float DX, float DY, float DZ, int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float* vx, float* vy, float* vz, float* sxx, float* syy, float* szz, float* sxy, float* syz, float* sxz, float* rjp, float* rkp, float* rip) {

  float dx = DT / DX;
  float dy = DT / DY;
  float dz = DT / DZ;

  int strip = nz2;
  int slice = nz2 * nx2;

  // note that I've changed (<=ny2) to (<ny2 - ny1) and so for nx2 and nz2
  for (int j = ny1; j < ny2 - ny1; j++) {
    for (int i = nx1; i < nx2 - nx1; i++) {
      for (int k = nz1; k < nz2 - nz1; k++) {
        int idx = j * slice + i * strip + k;


        float sxx_x = dx * (sxx[idx + strip] - sxx[idx]);
        float sxy_y = dy * (sxy[idx] - sxy[idx - slice]);
        float sxz_z = dz * (sxz[idx] - sxz[idx - 1]); /* R�ckw�rtsoperator */

        /* updating components of particle velocities */
        vx[idx] += ((sxx_x + sxy_y + sxz_z) / rip[idx]);

        float syy_y = dy * (syy[idx + slice] - syy[idx]);
        float sxy_x = dx * (sxy[idx] - sxy[idx - strip]);
        float syz_z = dz * (syz[idx] - syz[idx - 1]);

        vy[idx] += ((syy_y + sxy_x + syz_z) / rjp[idx]);

        float szz_z = dz * (szz[idx + 1] - szz[idx]);
        float sxz_x = dx * (sxz[idx] - sxz[idx - strip]);
        float syz_y = dy * (syz[idx] - syz[idx - slice]);

        vz[idx] += ((szz_z + sxz_x + syz_y) / rkp[idx]);
      }
    }
  }
}

float *read_float3(int n1, int n2, int n3, sf_file F) {
  assert(n1 == sf_n(sf_iaxa(F, 1)));
  assert(n2 == sf_n(sf_iaxa(F, 2)));
  assert(n3 == sf_n(sf_iaxa(F, 3)));

  float *p = sf_floatalloc(n1 * n2 * n3);
  sf_floatread(p, n1 * n2 * n3, F);

  return p;
}

void write_float3(int n1, int n2, int n3, float *p, const char *fn) {
  sf_axis ax1 = sf_maxa(n1, 0, 1);
  sf_axis ax2 = sf_maxa(n2, 0, 1);
  sf_axis ax3 = sf_maxa(n3, 0, 1);
  sf_file F = sf_output(fn);

  sf_oaxa(F, ax1, 1);
  sf_oaxa(F, ax2, 2);
  sf_oaxa(F, ax3, 3);

  sf_floatwrite(p, n1 * n2 * n3, F);

  sf_maxa_free(ax1);
  sf_maxa_free(ax2);
  sf_maxa_free(ax3);
}


int main(int argc, char **argv)
{
  sf_init(argc, argv);

  float DT, DX, DY, DZ;
  int NX, NY, NZ;

  assert(sf_getfloat("DT", &DT));
  assert(sf_getfloat("DX", &DX));
  assert(sf_getfloat("DY", &DY));
  assert(sf_getfloat("DZ", &DZ));
  assert(sf_getint("NX", &NX));
  assert(sf_getint("NY", &NY));
  assert(sf_getint("NZ", &NZ));


  sf_file f_vx = sf_input("vx");
  sf_file f_vy = sf_input("vy");
  sf_file f_vz = sf_input("vz");
  sf_file f_sxx = sf_input("sxx");
  sf_file f_syy = sf_input("syy");
  sf_file f_szz = sf_input("szz");
  sf_file f_sxy = sf_input("sxy");
  sf_file f_syz = sf_input("syz");
  sf_file f_sxz = sf_input("sxz");
  sf_file f_rip = sf_input("rip");
  sf_file f_rjp = sf_input("rjp");
  sf_file f_rkp = sf_input("rkp");

  float *vx = read_float3(NZ, NX, NY, f_vx);
  float *vy = read_float3(NZ, NX, NY, f_vy);
  float *vz = read_float3(NZ, NX, NY, f_vz);
  float *sxx = read_float3(NZ, NX, NY, f_sxx);
  float *syy = read_float3(NZ, NX, NY, f_syy);
  float *szz = read_float3(NZ, NX, NY, f_szz);
  float *sxy = read_float3(NZ, NX, NY, f_sxy);
  float *syz = read_float3(NZ, NX, NY, f_syz);
  float *sxz = read_float3(NZ, NX, NY, f_sxz);
  float *rip = read_float3(NZ, NX, NY, f_rip);
  float *rjp = read_float3(NZ, NX, NY, f_rjp);
  float *rkp = read_float3(NZ, NX, NY, f_rkp);

  fd2_addr1d(DT, DX, DY, DZ, 1, NX, 1, NY, 1, NZ, vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, rjp, rkp, rip);

  /// write to output file
  write_float3(NZ, NX, NY, vx, "vx1");
  write_float3(NZ, NX, NY, vy, "vy1");
  write_float3(NZ, NX, NY, vz, "vz1");

  free(vx);
  free(vy);
  free(vz);
  free(sxx);
  free(syy);
  free(szz);
  free(sxy);
  free(syz);
  free(sxz);
  free(rip);
  free(rjp);
  free(rkp);

  sf_close();
  printf("exit normally\n");

  return 0;
}
