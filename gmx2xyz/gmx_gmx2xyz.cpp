#include "gmx_gmx2xyz.hpp"

typedef struct {
    int bond[5];
    int nb;
} t_bonded;

typedef struct {
    std::string atomname;
    std::string resname;
    int  atomid;
} amoeba_parm;

void read_topol( const char *top_file, t_bonded *bonds) {
    char line[1024];
    gmx_bool has_bonding = FALSE;
    gmx_bool in_bonds = FALSE;
    
    FILE *f = fopen(top_file,"r");
    while (fgets(line,1024,f))
    {
        char char1[250];
        char char2[250];
        char junk[250];
        memset(char1,0,sizeof(char1));
        memset(char2,0,sizeof(char2));
        sscanf(line,"%s %s",char1,char2);
        
        if (in_bonds) {
            int atom1,atom2;
            atom1=-1;
            if (strncmp(char1,";",1) != 0 && strncmp(char1,"[",1) != 0)
            {
                atom1 = atoi(char1);
                atom2 = atoi(char2);
                bonds[atom1].bond[bonds[atom1].nb]=atom2;
                bonds[atom2].bond[bonds[atom2].nb]=atom1;
                bonds[atom1].nb++;
                bonds[atom2].nb++;
            }
        }
        if (strncmp(char1,"[",1) == 0) {
            if (strncmp(char2,"bond",4) == 0) {
                has_bonding = TRUE;
                in_bonds = TRUE;
            }
            else {in_bonds = FALSE;}
        }
    }
    if ( f == NULL ) {
        fprintf(stderr,"\nError opening %s",top_file);
        exit(1);
    }
    if ( !has_bonding) {
        fprintf(stderr,"\nError!  There are no bonding parameters in your topology file!  You may need to catenate your itp file into the topology.\n");
        exit(1);
    }
    
    fclose(f);
    return;
}

void ambernames(amoeba_parm &parm, std::string residue, std::string atomname, int atomid, std::unordered_map<std::string, std::string> &residue_name_map){
    if (atomname == "HN") {
        atomname = "H";
    }
    parm.atomname = atomname;
    parm.resname = residue_name_map[residue];
    parm.atomid = atomid;
    
    // SOL and GNP have naming issues
    if (parm.resname == "SOL") {
        if (atomname == "O") {
            parm.atomname = "OW";
        }
        else {
            parm.atomname = "HW";
        }
    }
    else if (parm.resname == "GNP") {
        if (strncmp(&atomname[3],"'",1) == 0) {
            std::stringstream rename;
            rename << &atomname[1] << atomname[0];
            parm.atomname = rename.str();
        }
    }
    // C-Terminal OC1 and OC2 are called OXT.  Need to relabel as OC
    if ( residue.find("C-Terminal") != std::string::npos) {
        if (parm.atomname == "OXT") {
            parm.atomname = "OC";
        }
    }
    // All the ion atomnames should not have a + or -.  Just remove last character
    if (parm.resname == "Na+" ) { parm.atomname="Na"; }
    if (parm.resname == "Cl-" ) { parm.atomname="Cl"; }
    if (parm.resname == "ZN" ) { parm.atomname="Zn"; }
    if (parm.resname == "CA" ) { parm.atomname="Ca"; }

    return;
}

void read_amoeba_parameters( const char *parm_name, std::unordered_map<std::string, std::unordered_map<std::string, amoeba_parm> > &mapping){
    std::unordered_map<std::string, std::string> residue_name_map =
    {
        // Unassigned are present in the version of amoeba.prm that I use but
        // I have no need of and/or don't know offhand the amber03 naming
        {"\"3'-Hydroxyl DNA        \"",""},
        {"\"3'-Hydroxyl RNA        \"",""},
        {"\"3'-Phosphate OP DNA    \"",""},
        {"\"3'-Phosphate OP RNA    \"",""},
        {"\"3'-Phosphate OS DNA    \"",""},
        {"\"3'-Phosphate OS RNA    \"",""},
        {"\"3'-Phosphate P DNA     \"",""},
        {"\"3'-Phosphate P RNA     \"",""},
        {"\"5'-Hydroxyl DNA        \"",""},
        {"\"5'-Hydroxyl RNA        \"",""},
        {"\"5'-Phosphate OP DNA    \"",""},
        {"\"5'-Phosphate OP RNA    \"",""},
        {"\"5'-Phosphate OS DNA    \"",""},
        {"\"5'-Phosphate OS RNA    \"",""},
        {"\"5'-Phosphate P DNA     \"",""},
        {"\"5'-Phosphate P RNA     \"",""},
        {"\"Acetyl N-Terminus      \"",""},
        {"\"Adenosine              \"",""},
        {"\"Alanine                \"","ALA"},
        {"\"Amide C-Terminus       \"",""},
        {"\"Arginine               \"","ARG"},
        {"\"Asparagine             \"","ASN"},
        {"\"Aspartic Acid          \"","ASP"},
        {"\"C-Terminal AIB         \"",""},
        {"\"C-Terminal ALA         \"","CALA"},
        {"\"C-Terminal ARG         \"","CARG"},
        {"\"C-Terminal ASN         \"","CASN"},
        {"\"C-Terminal ASP         \"","CASP"},
        {"\"C-Terminal CYS (SH)    \"","CCYS"},
        {"\"C-Terminal CYS (SS)    \"","CCYX"}, // There is no CCYM and CCYS has an extra H
        {"\"C-Terminal GLN         \"","CGLN"},
        {"\"C-Terminal GLU         \"","CGLU"},
        {"\"C-Terminal GLY         \"","CGLY"},
        {"\"C-Terminal HIS (+)     \"","CHIP"},
        {"\"C-Terminal HIS (HD)    \"","CHID"},
        {"\"C-Terminal HIS (HE)    \"","CHIE"},
        {"\"C-Terminal ILE         \"","CILE"},
        {"\"C-Terminal LEU         \"","CLEU"},
        {"\"C-Terminal LYS         \"","CLYS"},
        {"\"C-Terminal MET         \"","CMET"},
        {"\"C-Terminal ORN         \"",""},
        {"\"C-Terminal PHE         \"","CPHE"},
        {"\"C-Terminal PRO         \"","CPRO"},
        {"\"C-Terminal SER         \"","CSER"},
        {"\"C-Terminal THR         \"","CTHR"},
        {"\"C-Terminal TRP         \"","CTRP"},
        {"\"C-Terminal TYR         \"","CTYR"},
        {"\"C-Terminal VAL         \"","CVAL"},
        {"\"CNC                    \"","CNC"},
        {"\"Calcium Ion            \"","CA"},
        {"\"Chloride Ion           \"","Cl-"},
        {"\"Cysteine (SH)          \"","CYS"},
        {"\"Cystine (SS)           \"","CYM"},
        {"\"Cytidine               \"",""},
        {"\"DCN                    \"","DCN"},
        {"\"Deoxyadenosine         \"",""},
        {"\"Deoxycytidine          \"",""},
        {"\"Deoxyguanosine         \"",""},
        {"\"Deoxythymidine         \"",""},
        {"\"Formyl N-Terminus      \"",""},
        {"\"GNP                    \"","GNP"},
        {"\"Glutamic Acid          \"","GLU"},
        {"\"Glutamine              \"","GLN"},
        {"\"Glycine                \"","GLY"},
        {"\"Guanosine              \"",""},
        {"\"Histidine (+)          \"","HIP"},
        {"\"Histidine (HD)         \"","HID"},
        {"\"Histidine (HE)         \"","HIE"},
        {"\"Isoleucine             \"","ILE"},
        {"\"Leucine                \"","LEU"},
        {"\"Lysine                 \"","LYS"},
        {"\"Magnesium Ion          \"","MG"},
        {"\"Methionine             \"","MET"},
        {"\"MethylAlanine (AIB)    \"",""},
        {"\"N-MeAmide C-Terminus   \"","NME"},
        {"\"N-Terminal AIB         \"",""},
        {"\"N-Terminal ALA         \"","NALA"},
        {"\"N-Terminal ARG         \"","NARG"},
        {"\"N-Terminal ASN         \"","NASN"},
        {"\"N-Terminal ASP         \"","NASP"},
        {"\"N-Terminal CYS (SH)    \"","NCYS"},
        {"\"N-Terminal CYS (SS)    \"","NCYX"}, // There is no NCYM and NCYS has an extra H
        {"\"N-Terminal GLN         \"","NGLN"},
        {"\"N-Terminal GLU         \"","NGLU"},
        {"\"N-Terminal GLY         \"","NGLY"},
        {"\"N-Terminal HIS (+)     \"","NHIP"},
        {"\"N-Terminal HIS (HD)    \"","NHID"},
        {"\"N-Terminal HIS (HE)    \"","NHIE"},
        {"\"N-Terminal ILE         \"","NILE"},
        {"\"N-Terminal LEU         \"","NLEU"},
        {"\"N-Terminal LYS         \"","NLYS"},
        {"\"N-Terminal MET         \"","NMET"},
        {"\"N-Terminal ORN         \"",""},
        {"\"N-Terminal PHE         \"","NPHE"},
        {"\"N-Terminal PRO         \"","NPRO"},
        {"\"N-Terminal SER         \"","NSER"},
        {"\"N-Terminal THR         \"","NTHR"},
        {"\"N-Terminal TRP         \"","NTRP"},
        {"\"N-Terminal TYR         \"","NTYR"},
        {"\"N-Terminal VAL         \"","NVAL"},
        {"\"Ornithine              \"",""},
        {"\"Phenylalanine          \"","PHE"},
        {"\"Phosphodiester DNA     \"",""},
        {"\"Phosphodiester RNA     \"",""},
        {"\"Potassium Ion          \"","K"},
        {"\"Proline                \"","PRO"},
        {"\"Pyroglutamic Acid      \"",""},
        {"\"Serine                 \"","SER"},
        {"\"Sodium Ion             \"","Na+"},
        {"\"Threonine              \"","THR"},
        {"\"Tryptophan             \"","TRP"},
        {"\"Tyrosine               \"","TYR"},
        {"\"Uridine                \"",""},
        {"\"Valine                 \"","VAL"},
        {"\"Water                  \"","SOL"},
        {"\"Zinc Ion (+2)          \"","ZN"},
        // Accounting for the change in spacing in amoebapro13
        {"\"Acetyl N-Terminus\""      ,""},
        {"\"Alanine\""                ,"ALA"},
        {"\"Amide C-Terminus\""       ,""},
        {"\"Arginine\""               ,"ARG"},
        {"\"Asparagine\""             ,"ASN"},
        {"\"Aspartic Acid (COOH)\""   ,""},
        {"\"Aspartic Acid\""          ,"ASP"},
        {"\"C-Terminal AIB\""         ,""},
        {"\"C-Terminal ALA\""         ,"CALA"},
        {"\"C-Terminal ARG\""         ,"CARG"},
        {"\"C-Terminal ASH (COOH)\""  ,""},
        {"\"C-Terminal ASN\""         ,"CASN"},
        {"\"C-Terminal ASP\""         ,"CASP"},
        {"\"C-Terminal CYD (S-)\""    ,""},
        {"\"C-Terminal CYS (SH)\""    ,"CCYS"},
        {"\"C-Terminal CYX (SS)\""    ,"CCYX"},
        {"\"C-Terminal GLH (COOH)\""  ,""},
        {"\"C-Terminal GLN\""         ,"CGLN"},
        {"\"C-Terminal GLU\""         ,"CGLU"},
        {"\"C-Terminal GLY\""         ,"CGLY"},
        {"\"C-Terminal HIS (+)\""     ,"CHIP"},
        {"\"C-Terminal HIS (HD)\""    ,"CHID"},
        {"\"C-Terminal HIS (HE)\""    ,"CHIE"},
        {"\"C-Terminal ILE\""         ,"CILE"},
        {"\"C-Terminal LEU\""         ,"CLUE"},
        {"\"C-Terminal LYD (NH2)\""   ,""},
        {"\"C-Terminal LYS\""         ,"CLYS"},
        {"\"C-Terminal MET\""         ,"CMET"},
        {"\"C-Terminal ORN\""         ,""},
        {"\"C-Terminal PHE\""         ,"CPHE"},
        {"\"C-Terminal PRO\""         ,"CPRO"},
        {"\"C-Terminal SER\""         ,"CSER"},
        {"\"C-Terminal THR\""         ,"CTHR"},
        {"\"C-Terminal TRP\""         ,"CTRP"},
        {"\"C-Terminal TYD (O-)\""    ,""},
        {"\"C-Terminal TYR\""         ,"CTYR"},
        {"\"C-Terminal VAL\""         ,"CVAL"},
        {"\"CNC\""                    ,"CNC"},
        {"\"Calcium Ion\""            ,"CA"},
        {"\"Chloride Ion\""           ,"Cl-"},
        {"\"Cysteine (S-)\""          ,""},
        {"\"Cysteine (SH)\""          ,"CYS"},
        {"\"Cystine (SS)\""           ,"CYX"},
        {"\"Deprotonated N-Terminus\"",""},
        {"\"EtSCN\""                  ,"ESCN"},
        {"\"Formyl N-Terminus\""      ,""},
        {"\"Glutamic Acid (COOH)\""   ,"GLH"},
        {"\"Glutamic Acid\""          ,"GLU"},
        {"\"Glutamine\""              ,"GLN"},
        {"\"Glycine\""                ,"GLY"},
        {"\"Histidine (+)\""          ,"HIP"},
        {"\"Histidine (HD)\""         ,"HID"},
        {"\"Histidine (HE)\""         ,"HIE"},
        {"\"Isoleucine\""             ,"ILE"},
        {"\"Leucine\""                ,"LEU"},
        {"\"Lysine (NH2)\""           ,"LYP"},
        {"\"Lysine\""                 ,"LYS"},
        {"\"Magnesium Ion\""          ,"MG"},
        {"\"MeSCN\""                  ,"MSCN"},
        {"\"Methionine\""             ,"MET"},
        {"\"MethylAlanine (AIB)\""    ,""},
        {"\"N-MeAmide C-Terminus\""   ,"NME"},
        {"\"N-Terminal AIB\""         ,""},
        {"\"N-Terminal ALA\""         ,"NALA"},
        {"\"N-Terminal ARG\""         ,"NARG"},
        {"\"N-Terminal ASH (COOH)\""  ,""},
        {"\"N-Terminal ASN\""         ,"NASN"},
        {"\"N-Terminal ASP\""         ,"NASP"},
        {"\"N-Terminal CYD (S-)\""    ,""},
        {"\"N-Terminal CYS (SH)\""    ,"NCYS"},
        {"\"N-Terminal CYX (SS)\""    ,"NCYX"},
        {"\"N-Terminal GLH (COOH)\""  ,""},
        {"\"N-Terminal GLN\""         ,"NGLN"},
        {"\"N-Terminal GLU\""         ,"NGLU"},
        {"\"N-Terminal GLY\""         ,"NGLY"},
        {"\"N-Terminal HIS (+)\""     ,"NHIP"},
        {"\"N-Terminal HIS (HD)\""    ,"NHID"},
        {"\"N-Terminal HIS (HE)\""    ,"NHIE"},
        {"\"N-Terminal ILE\""         ,"NILE"},
        {"\"N-Terminal LEU\""         ,"NLEU"},
        {"\"N-Terminal LYD (NH2)\""   ,""},
        {"\"N-Terminal LYS\""         ,"NLYS"},
        {"\"N-Terminal MET\""         ,"NMET"},
        {"\"N-Terminal ORN\""         ,""},
        {"\"N-Terminal PHE\""         ,"NPHE"},
        {"\"N-Terminal PRO\""         ,"NPRO"},
        {"\"N-Terminal SER\""         ,"NSER"},
        {"\"N-Terminal THR\""         ,"NTHR"},
        {"\"N-Terminal TRP\""         ,"NTRP"},
        {"\"N-Terminal TYD (O-)\""    ,""},
        {"\"N-Terminal TYR\""         ,"NTYR"},
        {"\"N-Terminal VAL\""         ,"NVAL"},
        {"\"Ornithine\""              ,""},
        {"\"Phenylalanine\""          ,"PHE"},
        {"\"Potassium Ion\""          ,"K"},
        {"\"Proline\""                ,"PRO"},
        {"\"Protonated C-Terminus\""  ,""},
        {"\"Pyroglutamic Acid\""      ,""},
        {"\"Serine\""                 ,"SER"},
        {"\"Sodium Ion\""             ,"Na+"},
        {"\"Threonine\""              ,"THR"},
        {"\"Tryptophan\""             ,"TRP"},
        {"\"Tyrosine (O-)\""          ,""},
        {"\"Tyrosine\""               ,"TYR"},
        {"\"Valine\""                 ,"VAL"},
        {"\"Water\""                  ,"SOL"}
    };
    int n=0;
    std::string line;
    std::ifstream file(parm_name);
    if (file.is_open()){
        while (file.good()){
            getline(file,line);
            if (not line.empty()){
                std::string b,buffer,atomname,residue, satomid;
                int atomid;
                b = line.substr(0,7);
                if (b == "biotype"){
                    size_t quotestart = line.find('"');
                    size_t quoteend = line.find('"',quotestart+1) + 1;
                    residue = line.substr(quotestart,quoteend-quotestart);
                    std::string beforeQuotes = line.substr(0,quotestart-1);
                    std::string afterQuotes = line.substr(quoteend+1,line.length());
                    
                    std::stringstream ss(beforeQuotes+afterQuotes);
                    std::vector<std::string> tokens;
                    
                    while (ss >> buffer) {
                        tokens.push_back(buffer);
                    }
                    if (tokens.size() < 4) {
                        std::cerr << "Error reading in the parameter file.  Expected 3 columns prior to the group name in biotype and 1 column following. Offending line:\n" << line << std::endl;
                        exit(1);
                    }
                    atomname = tokens[2];
                    std::stringstream d3(tokens[3]);
                    d3 >> atomid;

                    amoeba_parm item;

                    ambernames(item,residue,atomname,atomid,residue_name_map);
                    mapping[item.resname].insert(std::make_pair(item.atomname,item));
                }
            }
        }
    }
    else if (!file){
        std::cerr << "\nError reading " << parm_name << std::endl;
        exit(1);
    }
    file.close();
    
    /*
    std::cout << "mymap contains:";
    for ( auto it = mapping.begin(); it != mapping.end(); ++it ) {
        std::cout << " " << it->first << " : ";
        for (auto jt = mapping[it->first].begin(); jt != mapping[it->first].end(); ++jt) {
            std::cout << jt->first << "(" << mapping[it->first][jt->first].atomname << ") ";
//            std::cout << mapping[it->first][jt->first].atomname << " ";
        }
    std::cout << std::endl;
    }
    */
    return;
};

int get_atomid( int &i, t_topology &top, std::unordered_map<std::string, std::unordered_map<std::string, amoeba_parm> > &atoms_types ) {
    int resid = top.atoms.atom[i].resind;
    std::string resname = *top.atoms.resinfo[resid].name;
    std::string atomname = *top.atoms.atomname[i];
    std::string original_atomname = *top.atoms.atomname[i];
    // If it's a C-Terminal or N-Terminal and backbone, use CXXX and NXXX respectively,
    // otherwise use XXX
    if ((resname.size() > 3) && (resname.substr(0,1) == "C" || resname.substr(0,1) == "N" )) {
        if (
            (atomname != "N"   && atomname != "CA"  && atomname != "C"   && atomname != "H"  &&
             atomname != "OC1" && atomname != "OC2" && atomname != "O"   && atomname != "HA" &&
             atomname != "H1"  && atomname != "H2"  && atomname != "H3"  )
            ) {
            resname = resname.substr(1,resname.size());
        }
    }
    while ( (atoms_types[resname][atomname].atomname != atomname) && atomname.size() > 1) {
//        std::cout << "\tChanging " << resname << " " << atomname << " to ";
        atomname = atomname.substr(0, atomname.size()-1);
//        std::cout << atomname << " (" << atomname.size() << ")" <<  std::endl;
    }
    if (atoms_types[resname][atomname].atomname == atomname) {
        if (atomname != original_atomname) {
//            std::cout << "Adding " << resname << " " << original_atomname << " from " << atomname << std::endl;
            // Add the shorted atomname so we don't have to go through this shortening again
            // next time
            atoms_types[resname].insert(std::make_pair(original_atomname,atoms_types[resname][atomname]));
        }
//        std::cout << "Assigned " << resname << " " << atomname << " to " << atoms_types[resname][atomname].atomid << std::endl;
        return atoms_types[resname][atomname].atomid;
    }
    else {
        std::cerr << "\nCannot find " << resname << " " << original_atomname << ", even after shortening to " << atomname << std::endl;
        // The problem is that terminals are missing side chains.  Either look at
        // non-terminal entries to get side chains, or find a way to add sidechains
        // to the terminal entries
        exit(1);
    }
    
    
}

void parm_order( t_topology &top, std::unordered_map<std::string, std::unordered_map<std::string, amoeba_parm> > &atoms_types, int &i_index, atom_id *index) {
    
    std::vector<int> orders;
    for (int n = 0; n<i_index; n++) {
        int i = index[n];
        int resid = top.atoms.atom[i].resind;
        std::string resname = *top.atoms.resinfo[resid].name;
        std::string atomname = *top.atoms.atomname[i];
        std::string original_atomname = *top.atoms.atomname[i];
        while ( (atoms_types[resname][atomname].atomname != atomname) && atomname.size() > 1) {
//            std::cout << "\tChanging " << atomname << " to ";
            atomname = atomname.substr(0, atomname.size()-1);
//            std::cout << atomname << " (" << atomname.size() << ")" << std::endl;
        }
        if (atoms_types[resname][atomname].atomname == atomname) {
            if (atomname != original_atomname) {
//                std::cout << "Adding " << resname << " " << original_atomname << " from " << atomname << std::endl;
                atoms_types[resname].insert(std::make_pair(original_atomname,atoms_types[resname][atomname]));
            }
//            std::cout << "Assigned " << resname << " " << atomname << " to " << atoms_types[resname][atomname].atomid << std::endl;
            orders.push_back(atoms_types[resname][atomname].atomid);
        }
        else {
            std::cerr << "Cannot find " << resname << " " << atomname << std::endl;
            // The problem is that terminals are missing side chains.  Either look at
            // non-terminal entries to get side chains, or find a way to add sidechains
            // to the terminal entries
            exit(1);
        }
    }
}

int gmx_gmx2xyz(int argc, char *argv[])
{
    const char      *desc[] = {
        "\tConvert from GROMACS format to Tinker XYZ format",
    };

    gmx_bool        bVerbose = FALSE;
    int             a1=0,a2=0;
    const char      *tpr_file, *top_file, *traj_file, *index_file, *out_file = NULL;
    char            tpr_title[256], xyz[256];
    int             i_index;
    atom_id         *ind_index;
    char            *gn_index;
    const char      *parm_name;
    t_bonded        *bonds;
    t_pargs         pa[] = {
        { "-v", FALSE, etBOOL,
            {&bVerbose}, "Be slightly more verbose"}
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        { efTOP, NULL, NULL, ffREAD },
        { efRND, "-a", "amoeba.prm", ffREAD },
        { efXYZ, "-o", NULL, ffWRITE },
        { efNDX, NULL, NULL, ffREAD }
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
    top_file = ftp2fn(efTOP, NFILE, fnm);
    init_top(&top);
    parm_name = opt2fn("-a", NFILE, fnm);
    out_file    = opt2fn("-o", NFILE, fnm);
    
    /* Open inputs */
    read_tps_conf(tpr_file, buffer, &top, &ePBC,
                  &xtop, NULL, box, TRUE);
    sfree(xtop);
    atoms = &top.atoms;
    
    /* Get index selection */
    printf("Select group for to save atom coordinates:\n");
    index_file = ftp2fn(efNDX, NFILE, fnm);
    get_index(atoms, index_file, 1, &i_index, &ind_index, &gn_index);
    snew(bonds,top.atoms.nr+1);
    read_topol( top_file, bonds);
    for (int i=0;i<top.atoms.nr;i++){
        if (strncmp(*top.atoms.atomname[i],"OW",2) == 0) {
            bonds[i+1].nb = 2;
            bonds[i+1].bond[0] = i+2;
            bonds[i+1].bond[1] = i+3;
        }
        else if (strncmp(*top.atoms.atomname[i],"HW1",3) ==0) {
            bonds[i+1].nb = 1;
            bonds[i+1].bond[0] = i-0;
        }
        else if (strncmp(*top.atoms.atomname[i],"HW2",3) ==0) {
            bonds[i+1].nb = 1;
            bonds[i+1].bond[0] = i-1;
        }
    }
    
    /* Read AMOEBA parameters */
    std::unordered_map<std::string, std::unordered_map<std::string, amoeba_parm> > atoms_types;
    read_amoeba_parameters(parm_name, atoms_types);
    /* Match AMBER parameters to AMOEBA parameters */
    std::vector<int> atomid_order(i_index,0);
    std::vector<std::string> atomname_order(i_index);
    for (int i=0;i<i_index;i++){
        atomid_order[i] = get_atomid(ind_index[i],top,atoms_types);
        // Also, since the names are the same every frame, just get those once
        std::stringstream name;
        name << *top.atoms.atomname[ind_index[i]] << "." << *top.atoms.resinfo[top.atoms.atom[ind_index[i]].resind].name;
        atomname_order[i] = name.str().c_str();
    }
    
    /* Read first frame */
    gmx_bool bHaveFirstFrame = read_first_frame(oenv, &status, traj_file, &fr, flags);
    if (bHaveFirstFrame) {
        set_trxframe_ePBC(&fr,ePBC);
    }
    
    /* read file and loop through frames */
    int frameN = 0;
    std::string outfile(out_file);
    do {
        std::stringstream outN;
        outN << outfile.substr(0,outfile.size()-4) << frameN << ".xyz";
        FILE * xyzout = ffopen(outN.str().c_str(),"w");
        fprintf(xyzout,"%7i  %s",i_index,outN.str().c_str());
        for (int n=0;n<i_index;n++){
            int i = ind_index[n];
            fprintf(xyzout, "\n%7i %9s %13.4f %13.4f %13.4f %7i ",
                    i+1,
                    atomname_order[n].c_str(),
                    fr.x[i][0]*10, fr.x[i][1]*10, fr.x[i][2]*10,
                    atomid_order[n]
                    );
            for (int j=0; j<bonds[i+1].nb; j++) {
                fprintf(xyzout,"%7i ",bonds[i+1].bond[j]);
            }
        }
        fclose(xyzout);
        frameN++;
    } while(read_next_frame(oenv, status, &fr));
    atomid_order.clear();
    atomname_order.clear();
    
    return 0;
}
