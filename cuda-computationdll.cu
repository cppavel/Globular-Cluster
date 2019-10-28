#include <iostream>
#include <algorithm>
#include "math.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "dll.h"

#define G  6.67408e-11

using namespace std;
namespace nbodydll
{
	__global__ void cuda_integrate(double* x, double* y, double* z, double* xacc, double* yacc, double* zacc, double* vx, double* vy, double*vz,
		double* md, double* m, double*datastorage, int cur_iteration, double timedelta, int cp)
	{
		int cur_body = blockIdx.x*blockDim.x + threadIdx.x;

		if (cur_body < cp)
		{
			xacc[cur_body] = 0;
			yacc[cur_body] = 0;
			zacc[cur_body] = 0;

			for (int i = 0; i < cp; i++)
			{
				if (i != cur_body)
				{
					xacc[cur_body] += (x[i] - x[cur_body]) * G * m[i] / ((md[i*cp + cur_body] * md[i*cp + cur_body]) * md[i*cp + cur_body]);
					yacc[cur_body] += (y[i] - y[cur_body]) * G* m[i] / ((md[i*cp + cur_body] * md[i*cp + cur_body]) * md[i*cp + cur_body]);
					zacc[cur_body] += (z[i] - z[cur_body]) * G * m[i] / ((md[i*cp + cur_body] * md[i*cp + cur_body]) * md[i*cp + cur_body]);
				}
			}

			x[cur_body] = x[cur_body] + vx[cur_body] * timedelta + xacc[cur_body] * timedelta*timedelta / 2;
			y[cur_body] = y[cur_body] + vy[cur_body] * timedelta + yacc[cur_body] * timedelta*timedelta / 2;
			z[cur_body] = z[cur_body] + vz[cur_body] * timedelta + zacc[cur_body] * timedelta*timedelta / 2;
			vx[cur_body] += xacc[cur_body] * timedelta;
			vy[cur_body] += yacc[cur_body] * timedelta;
			vz[cur_body] += zacc[cur_body] * timedelta;

			datastorage[(cur_iteration*cp + cur_body) * 3] = x[cur_body];
			datastorage[(cur_iteration*cp + cur_body) * 3 + 1] = y[cur_body];
			datastorage[(cur_iteration*cp + cur_body) * 3 + 2] = z[cur_body];
		}
	}

	__global__ void cuda_matrix_calculate(double* md, double* x, double* y, double* z, int cp)
	{
		int cur_body = blockIdx.x * blockDim.x + threadIdx.x;
		if (cur_body < cp)
		{
			for (int i = 0; i < cp; i++)
			{
				md[cur_body*cp + i] = md[i*cp + cur_body] = sqrt((x[cur_body] - x[i]) * (x[cur_body] - x[i]) + (y[cur_body] - y[i]) * (y[cur_body] - y[i])
					+ (z[cur_body] - z[i]) * (z[cur_body] - z[i]));
			}
		}

	}

	void cudaIntegrate(int* cp, double* x, double* y, double* z, double* vx, double* vy, double*vz, double* m, int iterations,
		double*datastorage, double* time, double startenergy, bool* merge, double* mergenergy)
	{

		//host
		double * md = new double[(*cp)*(*cp)];

		//device
		double * xdev;
		double *ydev;
		double *zdev;
		double *vxdev;
		double* vydev;
		double* vzdev;
		double* xaccdev;
		double* yaccdev;
		double* zaccdev;
		double* mddev;
		double* datastoragedev;
		double* mdev;

		cudaMalloc((double**)&xdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&ydev, (*cp) * sizeof(double));
		cudaMalloc((double**)&zdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&vxdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&vydev, (*cp) * sizeof(double));
		cudaMalloc((double**)&vzdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&xaccdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&yaccdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&zaccdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&mdev, (*cp) * sizeof(double));
		cudaMalloc((double**)&mddev, (*cp) * (*cp) * sizeof(double));
		cudaMalloc((double**)&datastoragedev, (*cp) * (iterations + 1) * 3 * sizeof(double));


		cudaMemcpy(xdev, x, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(ydev, y, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(zdev, z, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vxdev, vx, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vydev, vy, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(vzdev, vz, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(mdev, m, (*cp) * sizeof(double), cudaMemcpyHostToDevice);

		double deltatime = 0;
		int numberofblocks = (*cp) / 128 + 1;
		time[0] = 0;
		cuda_matrix_calculate << <numberofblocks, 128 >> > (mddev, xdev, ydev, zdev, (*cp));
		for (int i = 0; i <= iterations; i++)
		{
			if (i % 25 == 0)
			{
				cudaMemcpy(md, mddev, (*cp) * (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(vx, vxdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(vy, vydev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
				cudaMemcpy(vz, vzdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);

				double kinetic_energy = 0;
				for (int i = 0; i < (*cp); i++)
				{
					kinetic_energy += m[i] * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]) / 2;
				}
				double potential_energy = 0;
				for (int i = 0; i < (*cp); i++)
				{
					for (int j = i + 1; j < (*cp); j++)
					{
						potential_energy -= G*m[i] * m[j] / md[i*(*cp) + j];
					}
				}
				double new_total_energy = kinetic_energy + potential_energy + *mergenergy;
				double delta_energy = startenergy - new_total_energy;
				double coeff = sqrt(1 + delta_energy / kinetic_energy);
				for (int i = 0; i < (*cp); i++)
				{
					vx[i] *= coeff;
					vy[i] *= coeff;
					vz[i] *= coeff;
				}

				double x_impulse = 0;
				double y_impulse = 0;
				double z_impulse = 0;
				double sum_mass = 0;

				for (int i = 0; i < (*cp); i++)
				{
					x_impulse += m[i] * vx[i];
					y_impulse += m[i] * vy[i];
					z_impulse += m[i] * vz[i];
					sum_mass += m[i];
				}

				x_impulse = x_impulse / sum_mass;
				y_impulse = y_impulse / sum_mass;
				z_impulse = z_impulse / sum_mass;

				for (int i = 0; i < (*cp); i++)
				{
					vx[i] -= x_impulse;
					vy[i] -= y_impulse;
					vz[i] -= z_impulse;
				}

				double mindist = 1e290;
				int indexi = 0;
				int indexj = 0;
				for (int i = 0; i < (*cp); i++)
				{
					for (int j = i + 1; j < (*cp); j++)
					{
						if (md[i*(*cp) + j] < mindist)
						{
							indexi = i;
							indexj = j;
							mindist = md[i*(*cp) + j];
						}
					}
				}
				double maxv = 0;
				if (mindist < 1e9)
				{

					cudaMemcpy(md, mddev, (*cp)*(*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(x, xdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(y, ydev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(z, zdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(vx, vxdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(vy, vydev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(vz, vzdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
					cudaMemcpy(datastorage, datastoragedev, (*cp)*(iterations + 1) * 3 * sizeof(double), cudaMemcpyDeviceToHost);

					*merge = true;
					*mergenergy = *mergenergy - G*m[indexi] * m[indexj] / mindist;
					vx[indexi] = (m[indexi] * vx[indexi] + m[indexj] * vx[indexj]) / (m[indexi] + m[indexj]);
					vy[indexi] = (m[indexi] * vy[indexi] + m[indexj] * vy[indexj]) / (m[indexi] + m[indexj]);
					vz[indexi] = (m[indexi] * vz[indexi] + m[indexj] * vz[indexj]) / (m[indexi] + m[indexj]);

					x[indexi] = (m[indexi] * x[indexi] + m[indexj] * x[indexj]) / (m[indexi] + m[indexj]);
					y[indexi] = (m[indexi] * y[indexi] + m[indexj] * y[indexj]) / (m[indexi] + m[indexj]);
					z[indexi] = (m[indexi] * z[indexi] + m[indexj] * z[indexj]) / (m[indexi] + m[indexj]);

					m[indexi] = m[indexi] + m[indexj];

					for (int i = indexj; i < (*cp) - 1; i++)
					{
						x[i] = x[i + 1];
						y[i] = y[i + 1];
						z[i] = z[i + 1];
						vx[i] = vx[i + 1];
						vy[i] = vy[i + 1];
						vz[i] = vz[i + 1];
					}
					(*cp) = (*cp) - 1;

					for (int k = i; k < iterations + 1; k++)
					{
						time[k + 1] = time[i];
					}

					for (int k = i; k < iterations + 1; k++)
					{
						for (int j = 0; j < (*cp) + 1; j++)
						{
							datastorage[3 * (k*(*cp) + j)] = x[j];
							datastorage[3 * (k*(*cp) + j) + 1] = y[j];
							datastorage[3 * (k*(*cp) + j) + 2] = z[j];
						}

					}
					break;
				}

				maxv = sqrt(max(vx[indexi] * vx[indexi] + vy[indexi] * vy[indexi] + vz[indexi] * vz[indexi], vx[indexj] * vx[indexj] + vy[indexj] * vy[indexj] + vz[indexj] * vz[indexj]));
				deltatime = 1e-2*mindist / maxv;

				cudaMemcpy(vxdev, vx, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(vydev, vy, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
				cudaMemcpy(vzdev, vz, (*cp) * sizeof(double), cudaMemcpyHostToDevice);
			}
			if (i != 0)
			{
				time[i] = time[i - 1] + deltatime;
			}
			cuda_integrate << <numberofblocks, 128 >> > (xdev, ydev, zdev, xaccdev, yaccdev, zaccdev, vxdev, vydev, vzdev, mddev, mdev, datastoragedev, i, deltatime, (*cp));
			cuda_matrix_calculate << <numberofblocks, 128 >> > (mddev, xdev, ydev, zdev, (*cp));
		}
		if (!*merge)
		{
			cudaMemcpy(md, mddev, (*cp)*(*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(x, xdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(y, ydev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(z, zdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(vx, vxdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(vy, vydev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(vz, vzdev, (*cp) * sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(datastorage, datastoragedev, (*cp)*(iterations + 1) * 3 * sizeof(double), cudaMemcpyDeviceToHost);
		}


		delete[] md;

		cudaFree(xdev);
		cudaFree(ydev);
		cudaFree(zdev);
		cudaFree(vxdev);
		cudaFree(vydev);
		cudaFree(vzdev);
		cudaFree(xaccdev);
		cudaFree(yaccdev);
		cudaFree(zaccdev);
		cudaFree(mdev);
		cudaFree(mddev);
		cudaFree(datastoragedev);


	}
}
