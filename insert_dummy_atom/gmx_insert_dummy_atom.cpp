#include "gmx_insert_dummy_atom.hpp"

// Copied from confio.c
static void write_hconf_box(FILE *out, int pr, matrix box)
{
    char format[100];
    int  l;
    
    if (pr < 5)
    {
        pr = 5;
    }
    l = pr+5;
    
    if (box[XX][YY] || box[XX][ZZ] || box[YY][XX] || box[YY][ZZ] ||
        box[ZZ][XX] || box[ZZ][YY])
    {
        sprintf(format, "%%%d.%df%%%d.%df%%%d.%df"
                "%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df%%%d.%df\n",
                l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr, l, pr);
        fprintf(out, format,
                box[XX][XX], box[YY][YY], box[ZZ][ZZ],
                box[XX][YY], box[XX][ZZ], box[YY][XX],
                box[YY][ZZ], box[ZZ][XX], box[ZZ][YY]);
    }
    else
    {
        sprintf(format, "%%%d.%df%%%d.%df%%%d.%df\n", l, pr, l, pr, l, pr);
        fprintf(out, format,
                box[XX][XX], box[YY][YY], box[ZZ][ZZ]);
    }
}


int gmx_insert_dummy_atom(int argc, char *argv[])
{
    const char      *desc[] = {
        "\tAdd a dummy atom between atoms -a1 and -a2",
    };

    gmx_bool        bVerbose = FALSE;
    int             a1=-1,a2=-1;
	int 			insert_at=-1;
    const char      *tpr_file, *traj_file, *out_file;
    t_pargs         pa[] = {
        { "-a1", TRUE, etINT,
            {&a1}, "Starting atom for bond vector--ie: CD in CNC"},
        { "-a2", TRUE, etINT,
            {&a2}, "Ending atom for bond vector--ie: NE in CNC"},
        { "-insert_at", TRUE, etINT,
            {&insert_at}, "Desired index for dummy atom. Default = 0"},
        { "-v", FALSE, etBOOL,
            {&bVerbose}, "Be slightly more verbose"}
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        {efTRO, "-o","tilt",ffWRITE},
    };
#define NFILE asize(fnm)
#define NPA asize(pa)
    output_env_t    oenv;
    int             ngrps, nrefgrps;
    t_topology      top;
    t_atoms         *atoms=NULL;
    t_trxframe      fr,frout;
    t_trxstatus     *status;
    rvec            *xtop;
    matrix          box;
    int             ePBC;
    int             flags=TRX_READ_X;
    char            buffer[1024];
    int             ftp;
    FILE            *out_gro = NULL;
    t_trxstatus     *trxout = NULL;

    CopyRight(stderr,argv[0]);
    parse_common_args(&argc, argv, 
                      PCA_CAN_BEGIN | PCA_CAN_END | PCA_CAN_VIEW | 
                      PCA_TIME_UNIT | PCA_BE_NICE | PCA_CAN_TIME,
                      NFILE, fnm, NPA, pa, asize(desc), desc,
                      0, NULL, &oenv);

    /* Get inputs */
    tpr_file    = ftp2fn(efTPS, NFILE, fnm);
    traj_file   = opt2fn( "-f", NFILE, fnm);
    out_file    = opt2fn("-o", NFILE, fnm);
    ftp = fn2ftp(out_file);
    
    std::cout << "\n\n" << out_file << " " << ftp << std::endl;
    
    /* Open inputs */
    read_tps_conf(tpr_file, buffer, &top, &ePBC,
                  &xtop, NULL, box, TRUE);
    sfree(xtop);
    atoms = &top.atoms;
    gmx_conect gc = NULL;
    gc = gmx_conect_generate(&top);
    
    /* Make sure -a1 and -a2 are included and increment by -1 to match internal numbering */
    if ( a1<1 || a2<1 || a1==a2 || a1>top.atoms.nr || a2>top.atoms.nr ) {
        gmx_fatal(FARGS, "\nAtom numbers -a1 and -a2 defining the bond vector must be specified and different\n");
    }
    a1--; a2--;

    /* Read first frame */
    gmx_bool bHaveFirstFrame = read_first_frame(oenv, &status, traj_file, &fr, flags);
    if (bHaveFirstFrame) {
        set_trxframe_ePBC(&fr,ePBC);
    }
    
    if (ftp == efGRO) {
        out_gro = ffopen(out_file,"w");
    }
    else if (ftp == efXTC) {
        trxout = open_trx(out_file,"w");
    }
    
    /* read file and loop through frames */
    int frameN = 0;
    do {
        if (ftp == efGRO) {
            fprintf(out_gro,"Dummy atom inserted into %s, FRAME %i\n",traj_file,frameN);
            fprintf(out_gro,"%i\n",top.atoms.nr+1);
            float CD[3] = { fr.x[a1][0], fr.x[a1][1], fr.x[a1][2] };
            float NE[3] = { fr.x[a2][0], fr.x[a2][1], fr.x[a2][2] };
            float MP[3] ;
            for (int i=0;i<3;i++) {
                MP[i] = (CD[i]+NE[i])*0.5;
            }
            // GRO format:
            // RESID, RESNAME, ATOM, INDEX, X, Y, Z, vX, vY, vZ
            //"%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
            int index = 1;
			if (insert_at <= 0) {
            	fprintf(out_gro,"%5d%-5s%5s%5d%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
                    0,"TCHG","TCHG",0,MP[0],MP[1],MP[2],0.0f,0.0f,0.0f);
					index++;
			}

            /* Loop over atoms */
			int i;
			int resid_offset = 0;
            for (i=0;i<top.atoms.nr;i++){
				if (insert_at == index) {
            		fprintf(out_gro,"%5d%-5s%5s%5d%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
                    	top.atoms.atom[i-1].resind+2,"TCHG","TCHG",index,MP[0],MP[1],MP[2],0.0f,0.0f,0.0f);
					index++;
					resid_offset++;
				}
                // Ignoring velocities since I'm using this with mdrun -rerun for
                // force calculations only, which don't care about velocities
                fprintf(out_gro,"%5d%-5s%5s%5d%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
                        top.atoms.atom[i].resind+1 + resid_offset,
                        *top.atoms.resinfo[top.atoms.atom[i].resind].name,
                        *top.atoms.atomname[i],
                        index,
                        fr.x[i][0], fr.x[i][1], fr.x[i][2],
                        0.0f,0.0f,0.0f);
                index++;
                if (index > 99999) {
                    index = 0;
                }
            }
            /* Get box information */
            write_hconf_box(out_gro,1,box);
        }
        else if (ftp == efXTC) {
            float CD[3] = { fr.x[a1][0], fr.x[a1][1], fr.x[a1][2] };
            float NE[3] = { fr.x[a2][0], fr.x[a2][1], fr.x[a2][2] };
            rvec MP ;
            for (int i=0;i<3;i++) {
                MP[i] = (CD[i]+NE[i])*0.5;
            }
            rvec * newX = new rvec [top.atoms.nr+1];
			int i = 0;
			int offset = 0;
			if (insert_at <= 0) {
            	for (int j=0;j<3;j++) {
                	newX[i][j] = MP[j];
            	}
				offset++;
			}
            for (i=0; i<top.atoms.nr; i++) {
				if (insert_at == i) {
            		for (int j=0;j<3;j++) {
                		newX[i][j] = MP[j];
					}
					offset++;
				}
                for (int j=0;j<3;j++) {
                    newX[i+offset][j] = fr.x[i][j];
                }
            }
            frout = fr;
            frout.x = newX;
            frout.natoms++;
            write_trxframe(trxout,&frout,gc);
            delete[] newX;
        }
        frameN++;
    } while(read_next_frame(oenv, status, &fr));
    
    if (trxout) {
        close_trx(trxout);
    }
    if (out_gro) {
        ffclose(out_gro);
    }
    
    return 0;
}
