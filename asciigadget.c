#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"

FILE* filep;

int NDIM = 3;

static int n_type[6];

//TODO review reqd functions
void snapgadget(char*);

void fill_write_buffer(enum iofields, int*, int, int, int*, float*, float*, float*, int*);
size_t my_fwrite(void*, size_t, size_t, FILE *);

int main(int argc, char *argv[])
{
  char *fin, *fout;

  if(argc > 1)
    fin = argv[1];
  else
    {
      printf("no input file given!\n");
      exit(0);
    }
  if(argc > 2)
    fout = argv[2];
  else
    {
      printf("no output file given!\n");
      exit(0);
    }
  filep = fopen(fin, "r");

  snapgadget(fout);
  return (0);
}

void snapgadget(char* fout)
{
  int allocated_mass = 0;
  int allocated_pos = 0;
  int allocated_type = 0;
  int allocated_vel = 0;
  int allocated_key = 0;
  
  int *type_buf, *key_buf;
  float *pos_buf, *vel_buf;
  float *mass_buf;
  
  int nbody;
  float pmass;
  int i;
  float tsnap = 0.0;
  fscanf(filep, "%i %f", &nbody, &pmass);
  pos_buf = (float *) malloc(nbody * NDIM * sizeof(float));
  allocated_pos = 1;
  vel_buf = (float *) malloc(nbody * NDIM * sizeof(float));
  allocated_vel = 1;
  mass_buf = (float *) malloc(nbody * sizeof(float));
  allocated_mass = 1;
  key_buf = (int *) malloc(nbody * sizeof(int));
  allocated_key = 1;
  for(i = 0; i < nbody; i++)
    {
      if(pmass > 0.0)
	{
	  fscanf(filep, "%f %f %f %f %f %f", &pos_buf[i*NDIM], &pos_buf[i*NDIM+1], &pos_buf[i*NDIM+2], &vel_buf[i*NDIM], &vel_buf[i*NDIM+1], &vel_buf[i*NDIM+2]);
	  mass_buf[i] = (float) pmass;
	}
      else
	{
	  fscanf(filep, "%f %f %f %f %f %f %f", &mass_buf[i], &pos_buf[i*NDIM], &pos_buf[i*NDIM+1], &pos_buf[i*NDIM+2], &vel_buf[i*NDIM], &vel_buf[i*NDIM+1], &vel_buf[i*NDIM+2]);
	}
      
      key_buf[i] = i;
    }
  
  fclose(filep);
  
  printf("Read file.\n");
  
  All.BufferSize = 10000;
  
  //add a buffer for gadget type (here DM only)
  
  type_buf = (int *) malloc(nbody * sizeof(int));
  allocated_type = 1;
  int loop;
  for(loop = 0; loop < nbody; loop++)
    {
      type_buf[loop] = 1;
    }
  
  //now have all data read into bufs, proceed to write a gadget binary
  
  All.BufferSize = 100;
  
  //prepare the gadget header struct
  
  int type, bytes_per_blockelement, npart, nextblock, typelist[6];
  int n_for_this_task, ntask, n, p, pc, offset = 0, task;
  int blockmaxlen, ntot_type[6], nn[6];
  enum iofields blocknr;
  int blksize;
  //MPI_Status status;
  FILE *fd = 0;
  
#define SKIP  {my_fwrite(&blksize,sizeof(int),1,fd);}
  
  /* determine particle numbers of each tiype in file */ //TODO for now just 1 particle type, type 1 for halo DM particles
  
  for(n = 0; n < 6; n++)
    {
      ntot_type[n] = 0;
    }
  ntot_type[1] = nbody;
  
  for(n = 0; n < 6; n++)
    {
      n_type[n] = 0;
    }
  n_type[1] = nbody;
  
  
  /* fill file header */
  
  for(n = 0; n < 6; n++)
    {
      header.npart[n] = ntot_type[n];
      //header.npartTotal[n] = (unsigned int) ntot_type_all[n]; //TODO support for multiple files?
      header.npartTotal[n] = (unsigned int) ntot_type[n]; //TODO this supports single file only
      //header.npartTotalHighWord[n] = (unsigned int) (ntot_type_all[n] >> 32);
      header.npartTotalHighWord[n] = 0;
    }
  
  for(n = 0; n < 6; n++)
    {
      All.MassTable[n] = 0.0;
      if(n == 1)
	{
	  int write_masses = 0;
	  int m;
	  for(m = 1; m < nbody; m++)
	    {
	      if(mass_buf[m] != mass_buf[m-1])
		{
		  write_masses = 1;
		}
	    }
	  if(write_masses == 0)
	    All.MassTable[n] = mass_buf[0];
	  else
	    All.MassTable[n] = 0;
	}
    }
  
  for(n = 0; n < 6; n++)
    header.mass[n] = All.MassTable[n];
  
  All.Time = tsnap;
  header.time = All.Time;
  
  All.ComovingIntegrationOn = 0; //ZENO doesn't do cosmological ICs, I think - support for comoving integration has been REMOVED!
  header.redshift = 0;
  
  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;
  
  //TODO recheck these parameters
  All.NumFilesPerSnapshot = 1;
  All.BoxSize = 0.0; //not needed for non-periodic BCs, right?
  All.Omega0 = 0;
  All.OmegaLambda = 0; //comoving is off
  All.HubbleParam = 1.0;
  
  header.num_files = All.NumFilesPerSnapshot;
  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;
  
  /* Write header data to file */
  
  fd = fopen(fout, "w");
  
  blksize = sizeof(header);
  SKIP;
  my_fwrite(&header, sizeof(header), 1, fd);
  SKIP;
  
  /* Write particle data to file. */
  
  for(blocknr = 0; blocknr < IO_NBLOCKS; blocknr++)
    {
      if(blockpresent(blocknr))
	{ 
          bytes_per_blockelement = get_bytes_per_blockelement(blocknr);
	  blockmaxlen = ((int) (All.BufferSize * 1024 * 1024)) / bytes_per_blockelement;
	  npart = get_particles_in_block(blocknr, &typelist[0]);
	  
	  if(npart > 0)
	    {
              blksize = npart * bytes_per_blockelement;
	      SKIP;
	      
	      for(type = 0; type < 6; type++)
		{
		  if(typelist[type])
		    {
                      n_for_this_task = n_type[type];
		      while(1)
			{
                          pc = n_for_this_task;
			  
			  if(pc > blockmaxlen)
			    pc = blockmaxlen;
			  
			  CommBuffer = (void*) malloc(blockmaxlen*10000);

			  fill_write_buffer(blocknr, &offset, pc, type, type_buf, mass_buf, pos_buf, vel_buf, key_buf);
			  
			  my_fwrite(CommBuffer, bytes_per_blockelement, pc, fd);
			  
			  free(CommBuffer);
			  
			  n_for_this_task -= pc;
			  break;
			}
		    }
		}
	      SKIP;
	    }
	}
    }
  
  fclose(fd);
  
  //TODO re-check that these free's make sense with allocate's above
  if(allocated_mass)
    free(mass_buf);
  if(allocated_pos)
    free(pos_buf);
  if(allocated_vel)
    free(vel_buf);
  if(allocated_type)
    free(type_buf);
  if(allocated_key)
    free(key_buf);
}

/*! This function tells whether or not a given block in the output file is
 *  present, depending on the type of simulation run and the compile-time
 *  options. If one wants to add a new output-block, this function should be
 *  augmented accordingly.
 */
int blockpresent(enum iofields blocknr)
{
  //TODO review makefile for these flags

#ifndef OUTPUTPOTENTIAL
  if(blocknr == IO_POT)
    return 0;
#endif

#ifndef OUTPUTACCELERATION
  if(blocknr == IO_ACCEL)
    return 0;
#endif

#ifndef OUTPUTCHANGEOFENTROPY
  if(blocknr == IO_DTENTR)
    return 0;
#endif

#ifndef OUTPUTTIMESTEP
  if(blocknr == IO_TSTP)
    return 0;
#endif

  return 1;			/* default: present */
}

/*! This function tells the size of one data entry in each of the blocks
 *  defined for the output file. If one wants to add a new output-block, this
 *  function should be augmented accordingly.
 */
int get_bytes_per_blockelement(enum iofields blocknr)
{
  int bytes_per_blockelement = 0;

  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
      bytes_per_blockelement = 3 * sizeof(float);
      break;

    case IO_ID:
#ifdef LONGIDS 
      bytes_per_blockelement = sizeof(long long);
#else
      bytes_per_blockelement = sizeof(int);
#endif
      break;

    case IO_MASS:
    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_POT:
    case IO_DTENTR:
    case IO_TSTP:
      bytes_per_blockelement = sizeof(float);
      break;
    }

  return bytes_per_blockelement;
}

/*! This function determines how many particles there are in a given block,
 *  based on the information in the header-structure.  It also flags particle
 *  types that are present in the block in the typelist array. If one wants to
 *  add a new output-block, this function should be augmented accordingly.
 */
int get_particles_in_block(enum iofields blocknr, int *typelist)
{
  int i, nall, ntot_withmasses, ngas, nstars;

  nall = 0;
  ntot_withmasses = 0;

  for(i = 0; i < 6; i++)
    {
      typelist[i] = 0;

      if(header.npart[i] > 0)
	{
	  nall += header.npart[i];
	  typelist[i] = 1;
	}

      if(All.MassTable[i] == 0)
	ntot_withmasses += header.npart[i];
    }

  ngas = header.npart[0];
  nstars = header.npart[4];


  switch (blocknr)
    {
    case IO_POS:
    case IO_VEL:
    case IO_ACCEL:
    case IO_TSTP:
    case IO_ID:
    case IO_POT:
      return nall;
      break;

    case IO_MASS:
      for(i = 0; i < 6; i++)
	{
	  typelist[i] = 0;
	  if(All.MassTable[i] == 0 && header.npart[i] > 0)
	    typelist[i] = 1;
	}
      return ntot_withmasses;
      break;

    case IO_U:
    case IO_RHO:
    case IO_HSML:
    case IO_DTENTR:
      for(i = 1; i < 6; i++)
	typelist[i] = 0;
      return ngas;
      break;
    }

  exit(212);
  return 0;
}

/*! This catches I/O errors occuring for my_fwrite(). In this case we
 *  better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured: %s\n", ThisTask, strerror(errno));
      fflush(stdout);
      exit(777);
    }
  return nwritten;
}

/*! This function fills the write buffer with particle data. New output blocks
 *  can in principle be added here.
 */
void fill_write_buffer(enum iofields blocknr, int *startindex, int pc, int type, int* type_buf, float* mass_buf, float* pos_buf, float* vel_buf, int* key_buf)
{
  
  int n, k, pindex;
  float *fp;

#ifdef LONGIDS
  long long *ip;
#else
  int *ip;
#endif

  double dt_gravkick, dt_hydrokick, a3inv = 1, fac1, fac2;

  a3inv = fac1 = fac2 = 1;

  fp = CommBuffer;
  ip = CommBuffer;

  pindex = 0;

  switch (blocknr)
    {
    case IO_POS:		/* positions */
      printf("writing pos\n");
      for(n = 0; n < pc; pindex++)
	{
	  //if(type_buf[pindex] == type)
	  {
	    //printf("  n pc pindex %i %i %i\n",n,pc,pindex);
	    for(k = 0; k < 3; k++)
	      {
		fp[k] = pos_buf[pindex * NDIM + k];
	      }
	    n++;
	    fp += 3;
	  }
	}
      printf("pos written\n");
      break;

    case IO_VEL:		/* velocities */
      printf("writing vel\n");
      for(n = 0; n < pc; pindex++)
	//if(type_buf[pindex] == type)
	{
	  //printf("n pc pindex %i %i %i\n",n,pc,pindex);
	  for(k = 0; k < 3; k++)
	    {
	      //printf("  %f\n", vel_buf[pindex*NDIM+k]);
	      fp[k] = vel_buf[pindex * NDIM + k];
	    }
	  
	  n++;
	  fp += 3;
	}
      printf("vel written\n");
      break;

    case IO_ID:		/* particle ID */
      printf("writing ID\n");
      for(n = 0; n < pc; pindex++)
	//if(type_buf[pindex] == type)
	  {
	    *ip++ = key_buf[pindex];
	    n++;
	  }
      printf("ID written\n");
      break;

    case IO_MASS:		/* particle mass */
      printf("writing mass\n");
      for(n = 0; n < pc; pindex++)
	//if(type_buf[pindex] == type)
	  {
	    *fp++ = mass_buf[pindex];
	    n++;
	  }
      printf("mass written\n");
      break;


    }

  *startindex = pindex;
}
