#include <rsf.h>
#include <assert.h>

/**
 * Fuse the arrays. Please mind the order of each array
 * (vx, vy, vz) => Fv
 *  0   1    2
 *
 * (sxx, syy, szz, sxy, syz, sxz) => Fs
 *   0    1    2    3    4    5
 *
 * (rip, rjp, rkp) => Frp
 *   0    1    2
 */
void fd2_fuse(float DT, float DX, float DY, float DZ, int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float *Fv, float *Fs, float *Frp) {
  float dx = DT / DX;
  float dy = DT / DY;
  float dz = DT / DZ;

  int strip = nz2;
  int slice = nz2 * nx2;

  // note that I've changed (<=ny2) to (<ny2 - ny1) and so for nx2 and nz2
  for (int j = ny1; j < ny2 - ny1; j++) {
    for (int i = nx1; i < nx2 - nx1; i++) {
      for (int k = nz1; k < nz2 - nz1; k++) {
        /// each fused array maintains its own index
        int idx = j * slice + i * strip + k;


        float sxx_x = dx * (Fs[(idx + strip) * 6 + 0] - Fs[idx * 6]);
        float sxy_y = dy * (Fs[idx * 6 + 3] - Fs[(idx - slice) * 6 + 3]);
        float sxz_z = dz * (Fs[idx * 6 + 5] - Fs[(idx - 1) * 6 + 5]); /* R�ckw�rtsoperator */

        /* updating components of particle velocities */
        Fv[idx * 3 + 0] += ((sxx_x + sxy_y + sxz_z) / Frp[idx * 3 + 0]);

        float syy_y = dy * (Fs[(idx + slice) * 6 + 1] - Fs[idx * 6 + 1]);
        float sxy_x = dx * (Fs[idx * 6 + 3] - Fs[(idx - strip) * 6 + 3]);
        float syz_z = dz * (Fs[idx * 6 + 4] - Fs[(idx - 1) * 6 + 4]);

        Fv[idx * 3 + 1] += ((syy_y + sxy_x + syz_z) / Frp[idx * 3 + 1]);

        float szz_z = dz * (Fs[(idx + 1) * 6 + 2] - Fs[idx * 6 + 2]);
        float sxz_x = dx * (Fs[idx * 6 + 5] - Fs[(idx - strip) * 6 + 5]);
        float syz_y = dy * (Fs[idx * 6 + 4] - Fs[(idx - slice) * 6 + 4]);

        Fv[idx * 3 + 2] += ((szz_z + sxz_x + syz_y) / Frp[idx * 3 + 2]);
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

float *fuse_3(const float *a, const float *b, const float *c, int size) {
  float *p = sf_floatalloc(size * 3);
  for (int i = 0; i < size; i++) {
    p[i * 3 + 0] = a[i];
    p[i * 3 + 1] = b[i];
    p[i * 3 + 2] = c[i];
  }

  return p;
}

float *fuse_6(const float *a0, const float *a1, const float *a2, const float *a3, const float *a4, const float *a5, int size) {
  float *p = sf_floatalloc(size * 6);
  for (int i = 0; i < size; i++) {
    p[i * 6 + 0] = a0[i];
    p[i * 6 + 1] = a1[i];
    p[i * 6 + 2] = a2[i];
    p[i * 6 + 3] = a3[i];
    p[i * 6 + 4] = a4[i];
    p[i * 6 + 5] = a5[i];
  }

  return p;
}


void unfuse_3(const float *p, float *a, float *b, float *c, int size) {
  for (int i = 0; i < size; i++) {
    a[i] = p[i * 3 + 0];
    b[i] = p[i * 3 + 1];
    c[i] = p[i * 3 + 2];
  }
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

  float *Fv = fuse_3(vx, vy, vz, NX * NY * NZ);
  float *Fs = fuse_6(sxx, syy, szz, sxy, syz, sxz, NX * NY * NZ);
  float *Frp = fuse_3(rip, rjp, rkp, NX * NY * NZ);

  fd2_fuse(DT, DX, DY, DZ, 1, NX, 1, NY, 1, NZ, Fv, Fs, Frp);

  /// write to output file
  unfuse_3(Fv, vx, vy, vz, NX * NY * NZ);
  write_float3(NZ, NX, NY, vx, "vxf");
  write_float3(NZ, NX, NY, vy, "vyf");
  write_float3(NZ, NX, NY, vz, "vzf");

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
  free(Fv);
  free(Fs);
  free(Frp);

  sf_close();
  printf("exit normally\n");

  return 0;
}
