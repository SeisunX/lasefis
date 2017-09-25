#include <rsf.h>
#include <assert.h>


/**
 * Address the array in 3D ways, which is a good way for demonstrating the correctness, but might not be a CPE-friendly way
 */
void fd2_addr3d(float DT, float DX, float DY, float DZ, int nx1, int nx2, int ny1, int ny2, int nz1, int nz2, float*** vx, float*** vy, float*** vz, float*** sxx, float*** syy, float*** szz, float*** sxy, float*** syz, float*** sxz, float*** rjp, float*** rkp, float*** rip) {

  float dx = DT / DX;
  float dy = DT / DY;
  float dz = DT / DZ;

  // note that I've changed (<=ny2) to (<ny2 - ny1) and so for nx2 and nz2
  for (int j = ny1; j < ny2 - ny1; j++) {
    for (int i = nx1; i < nx2 - nx1; i++) {
      for (int k = nz1; k < nz2 - nz1; k++) {

        float sxx_x = dx * (sxx[j][i + 1][k] - sxx[j][i][k]);
        float sxy_y = dy * (sxy[j][i][k] - sxy[j - 1][i][k]);
        float sxz_z = dz * (sxz[j][i][k] - sxz[j][i][k - 1]); /* R�ckw�rtsoperator */

        /* updating components of particle velocities */
        vx[j][i][k] += ((sxx_x + sxy_y + sxz_z) / rip[j][i][k]);

        float syy_y = dy * (syy[j + 1][i][k] - syy[j][i][k]);
        float sxy_x = dx * (sxy[j][i][k] - sxy[j][i - 1][k]);
        float syz_z = dz * (syz[j][i][k] - syz[j][i][k - 1]);

        vy[j][i][k] += ((syy_y + sxy_x + syz_z) / rjp[j][i][k]);

        float szz_z = dz * (szz[j][i][k + 1] - szz[j][i][k]);
        float sxz_x = dx * (sxz[j][i][k] - sxz[j][i - 1][k]);
        float syz_y = dy * (syz[j][i][k] - syz[j - 1][i][k]);

        vz[j][i][k] += ((szz_z + sxz_x + syz_y) / rkp[j][i][k]);
      }
    }
  }
}

float ***read_float3(int n1, int n2, int n3, sf_file F) {
  assert(n1 == sf_n(sf_iaxa(F, 1)));
  assert(n2 == sf_n(sf_iaxa(F, 2)));
  assert(n3 == sf_n(sf_iaxa(F, 3)));

  float ***p = sf_floatalloc3(n1, n2, n3);
  sf_floatread(p[0][0], n1 * n2 * n3, F);

  return p;
}

void write_float3(int n1, int n2, int n3, float ***p, const char *fn) {
  sf_axis ax1 = sf_maxa(n1, 0, 1);
  sf_axis ax2 = sf_maxa(n2, 0, 1);
  sf_axis ax3 = sf_maxa(n3, 0, 1);
  sf_file F = sf_output(fn);

  sf_oaxa(F, ax1, 1);
  sf_oaxa(F, ax2, 2);
  sf_oaxa(F, ax3, 3);

  sf_floatwrite(p[0][0], n1 * n2 * n3, F);

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

  float ***vx = read_float3(NZ, NX, NY, f_vx);
  float ***vy = read_float3(NZ, NX, NY, f_vy);
  float ***vz = read_float3(NZ, NX, NY, f_vz);
  float ***sxx = read_float3(NZ, NX, NY, f_sxx);
  float ***syy = read_float3(NZ, NX, NY, f_syy);
  float ***szz = read_float3(NZ, NX, NY, f_szz);
  float ***sxy = read_float3(NZ, NX, NY, f_sxy);
  float ***syz = read_float3(NZ, NX, NY, f_syz);
  float ***sxz = read_float3(NZ, NX, NY, f_sxz);
  float ***rip = read_float3(NZ, NX, NY, f_rip);
  float ***rjp = read_float3(NZ, NX, NY, f_rjp);
  float ***rkp = read_float3(NZ, NX, NY, f_rkp);

  fd2_addr3d(DT, DX, DY, DZ, 1, NX, 1, NY, 1, NZ, vx, vy, vz, sxx, syy, szz, sxy, syz, sxz, rjp, rkp, rip);

  /// write to output file
  write_float3(NZ, NX, NY, vx, "vx3");
  write_float3(NZ, NX, NY, vy, "vy3");
  write_float3(NZ, NX, NY, vz, "vz3");

  free(**vx);  free(*vx);  free(vx);
  free(**vy);  free(*vy);  free(vy);
  free(**vz);  free(*vz);  free(vz);
  free(**sxx); free(*sxx); free(sxx);
  free(**syy); free(*syy); free(syy);
  free(**szz); free(*szz); free(szz);
  free(**sxy); free(*sxy); free(sxy);
  free(**syz); free(*syz); free(syz);
  free(**sxz); free(*sxz); free(sxz);
  free(**rip); free(*rip); free(rip);
  free(**rjp); free(*rjp); free(rjp);
  free(**rkp); free(*rkp); free(rkp);

  sf_close();
  printf("exit normally\n");

  return 0;
}
