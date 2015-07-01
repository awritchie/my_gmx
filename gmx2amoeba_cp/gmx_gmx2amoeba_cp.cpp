#include "gmx_gmx2amoeba_cp.hpp"

typedef struct {
    int bond[5];
    int nb;
    int near[10];
    int ib;
} t_bonded;

typedef struct {
    int atomid;
    int atomicnum;
    int nval;
    float alpha1;
    float beta1;
    float alpha2;
    float beta2;
    float alpha3;
    float beta3;
} valence;

typedef struct {
    int atomid;
    std::vector<int> bonded_id;
    int nbonded_ids;
    double q;
} qq;

typedef struct {
    std::string atomname;
    std::string resname;
    int  atomid;
    double q;
    std::vector<qq> qs;
    valence val;
} amoeba_parm;

void fill_val(int atomicnum, int nval, float alpha1, float beta1, float alpha2, float beta2, valence &item) {
    item.atomicnum = atomicnum;
    item.nval = nval;
    item.alpha1 = alpha1;
    item.beta1 = beta1;
    item.alpha2 = alpha2;
    item.beta2 = beta2;
    return;
}

void chgpen_field( const double *xyz, const double &r, const double &q, const double &z, const double &alpha, double *chgpen_field0){
    double rr = 1. / r;
    double rr2 = rr * rr;
    double rr3 = rr2 * rr;
    double tmp_exp = (alpha == 0) ? 0 : exp(-alpha*r);
    chgpen_field0[0];
    for (int j=0; j<3; j++){
        chgpen_field0[j] =    q * xyz[j] * rr3
                            - q * xyz[j] * rr3 * tmp_exp
                            - q * xyz[j] * rr2 * alpha * tmp_exp
                            + z * xyz[j] * rr3 * tmp_exp
                            + z * xyz[j] * rr2 * alpha * tmp_exp;
        chgpen_field0[j] *= cfac;
    }
    return;
}
double chgpen_switch_function( float r, float r_lo, float r_hi ){
    if (r <= r_lo){ return 1; }
    else if (r >= r_hi){ return 0;}
    else {
        float ratio = (r_hi - r)/(r_hi - r_lo);
        return 10*pow(ratio,3) - 15*pow(ratio,4) + 6*pow(ratio,5);
    }
}
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

void str2val( std::string &a, int &b) {
    std::stringstream a1(a);
    a1 >> b;
}
void str2val( std::string &a, float &b) {
    std::stringstream a1(a);
    a1 >> b;
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
        {"\"Zinc Ion (+2)          \"","ZN"}
    };
    // First, let's fine the charge on each atomid
    std::vector<float> charges(5000,0);
    std::vector<std::vector<qq> > id_charges(5000);
    std::string line;
    std::ifstream ffile(parm_name);
    if (ffile.is_open()){
        while (ffile.good()){
            getline(ffile,line);
            if (not line.empty()){
                std::string mp = line.substr(0,9);
                if ( mp == "multipole" ) {
                    double test = 0;
                    double ia=test, ib=test, ic=test, id=test, iq=test;
                    std::stringstream lline(line);
                    lline >> mp >> ia >> ib >> ic >> id >> iq;
//                    std::cout << line << std::endl;

                    if ( iq != test ) {
//                        std::cout << ia << " " << ib << " " << ic << " " << id << " " << iq << std::endl;
                    }
                    else if ( id != test ) {
                        iq = id;
//                        std::cout << ia << " " << ib << " " << ic << " " << iq << std::endl;
                    }
                    else if ( ic != test ) {
                        iq = ic;
//                        std::cout << ia << " " << ib << " " << iq << std::endl;
                    }
                    else if ( ib != test ) {
                        iq = ib;
//                        std::cout << ia << " " << iq << std::endl;
                    }
                    ia = abs(ia);
                    ib = abs(ib);
                    ic = abs(ic);
                    id = abs(id);
                    while ( (int)charges.size() < ia ) {
                        charges.push_back(0);
                    }
                    charges[ia] = iq;

                    while ( (int)id_charges.size() < ia) {
                        std::vector<qq> tmp;
                        id_charges.push_back(tmp);
                    }
                    qq multipole;
                    multipole.bonded_id = std::vector<int> (3,0);
                    std::stringstream llline(line);
                    llline >> mp >> ia;
                    std::vector<double> values;
                    double value;
                    while (llline >> value){
                        values.push_back(value);
                    }
                    for (int j=0;j<(int)values.size()-1;j++){
                        multipole.bonded_id[j]=values[j];
                    }
                    multipole.atomid = ia;
                    multipole.q = values[values.size()-1]; // charge is the last entry always
                    multipole.nbonded_ids = (int)multipole.bonded_id.size();
                    id_charges[ia].push_back(multipole);
                }
            }
        }
    }
    else if (!ffile){
        std::cerr << "\nError reading " << parm_name << std::endl;
        exit(1);
    }
    ffile.close();
    
    // Now, let's find the number of valence electrons on each atom
    std::vector<valence> from_table(36);
    // fill_val(int atomicnum, int nval, float alpha1, float beta1, float alpha2, float beta2, &valence item)
    for (int i=0;i<(int)from_table.size();i++){
        fill_val(i,0,0,0,0,0,from_table[i]);
    }
    fill_val(1,  1, 2., (1.999+2.010+2.004)/3., 3.3, (2.893,3.041,3.268)/3.,from_table[1]);  // H
    fill_val(6,  4, 4., (2.646+2.708+2.685)/3., 3.8, (3.021+2.771+2.724)/3.,from_table[6]);  // C
    fill_val(7,  5, 5., (3.097+3.072+3.054)/3., 3.1, (2.767+2.759+2.733)/3.,from_table[7]);  // N
    fill_val(8,  6, 6., (3.661+4.282+4.469)/3., 3.5, (3.083+3.350+3.170)/3.,from_table[8]);  // O
    fill_val(15, 5, 5.,                  2.360,  .0,                     .0,from_table[15]); // P
    fill_val(16, 6, 6.,       (2.770+2.381)/2.,  .0,                     .0,from_table[16]); // S
    fill_val(9,  7, 7.,                  4.275,  .0,                     .0,from_table[9]);  // F
    fill_val(17, 7, 7.,                  2.830,  .0,                     .0,from_table[17]); // Cl
    fill_val(35, 7, 7.,                  2.564,  .0,                     .0,from_table[35]); // Br

    std::vector<valence> val_parm(5000);
    std::ifstream fffile(parm_name);
    if (fffile.is_open()){
        while (fffile.good()){
            getline(fffile,line);
            if (not line.empty()){
                std::string mp = line.substr(0,4);
                if ( mp == "atom" ) {
                    int atomid, atomicn;
                    std::vector<std::string> sepline;
                    std::string item;
                    std::stringstream ll(line);
                    while ( ll >> item ) {
                        sepline.push_back(item);
                    }
                    str2val(sepline[1], atomid);
                    str2val(sepline[(int)sepline.size() -3], atomicn);
    
                    int nval;
                    float alpha1, beta1, alpha2, beta2;
                    if (atomicn <= 2){
                        nval = atomicn;
                    }
                    else if (atomicn <= 10){
                        nval = atomicn - 2;
                    }
                    else if (atomicn <= 18){
                        nval = atomicn - 10;
                    }
                    else if (atomicn <= 20){
                        nval = atomicn - 18;
                    }
                    else if (atomicn <= 36){
                        nval = max(2,atomicn - 28);
                    }
                    else if (atomicn <= 38){
                        nval = atomicn - 36;
                    }
                    else if (atomicn <= 54){
                        nval = max(2,atomicn - 48);
                    }
                    else if (atomicn <=56){
                        nval = atomicn - 54;
                    }
                    
                    if (sepline[3].compare(sepline[3].size()-1,1,"+") == 0){
                        nval--;
                    }
                    else if (sepline[3].compare(sepline[3].size()-1,1,"-") == 0){
                        nval++;
                    }
                    alpha1 = max(2,nval);
                    
                    val_parm[atomid].atomid = atomid;
                    val_parm[atomid].atomicnum = atomicn;
                    val_parm[atomid].nval = nval;
                    val_parm[atomid].alpha1 = alpha1;
                    val_parm[atomid].beta1 = from_table[atomicn].beta1;
                    val_parm[atomid].alpha2 = from_table[atomicn].alpha2;
                    val_parm[atomid].beta2 = from_table[atomicn].beta2;
                    val_parm[atomid].alpha3 = (from_table[atomicn].alpha2 == 0) ? alpha1 : from_table[atomicn].alpha2;
                    val_parm[atomid].beta3  = (from_table[atomicn].beta2  == 0) ? alpha1 : from_table[atomicn].beta2;

                    //std::cout << "chgpenprm " << atomid << " " << nval << " " << val_parm[atomid].alpha1 << " " << val_parm[atomid].beta1 << std::endl;
                    //std::cout << nval << " " << val_parm[atomid].alpha1 << " " << val_parm[atomid].alpha2 << " " << val_parm[atomid].alpha3 << std::endl;
                }
            }
        }
    }
    else if (!fffile){
        std::cerr << "\nError reading " << parm_name << std::endl;
        exit(1);
    }
    fffile.close();
    

    int n=0;
    std::ifstream file(parm_name);
    if (file.is_open()){
        while (file.good()){
            getline(file,line);
            if (not line.empty()){
                std::string b,buffer,atomname,residue, satomid;
                int atomid;
                b = line.substr(0,7);
                if (b == "biotype"){
                    atomname = line.substr(15,5);
                    std::stringstream d1(atomname);
                    d1 >> atomname;
                    residue = line.substr(22,25);
                    // I want to keep the whitespace here only
                    satomid = line.substr(50,4);
                    std::stringstream d3(satomid);
                    d3 >> atomid;
                    amoeba_parm item;
                    ambernames(item,residue,atomname,atomid,residue_name_map);
                    mapping[item.resname].insert(std::make_pair(item.atomname,item));
                    mapping[item.resname][item.atomname].q = charges[atomid];
                    mapping[item.resname][item.atomname].qs = id_charges[atomid];
                    mapping[item.resname][item.atomname].val = val_parm[atomid];
                    //std::cout << residue << " "<< atomname << " "<< atomid << " "<< item.resname << " " << item.atomname << std::endl;
                }
            }
        }
    }
    else if (!file){
        std::cerr << "\nError re-reading " << parm_name << std::endl;
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
    }*/
    return;
};

amoeba_parm get_atomid( int &i, t_topology &top, std::unordered_map<std::string, std::unordered_map<std::string, amoeba_parm> > &atoms_types ) {
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
        return atoms_types[resname][atomname];
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
            std::cout << "\tChanging " << atomname << " to ";
            atomname = atomname.substr(0, atomname.size()-1);
            std::cout << atomname << " (" << atomname.size() << ")" << std::endl;
        }
        if (atoms_types[resname][atomname].atomname == atomname) {
            if (atomname != original_atomname) {
                std::cout << "Adding " << resname << " " << original_atomname << " from " << atomname << std::endl;
                atoms_types[resname].insert(std::make_pair(original_atomname,atoms_types[resname][atomname]));
            }
            std::cout << "Assigned " << resname << " " << atomname << " to " << atoms_types[resname][atomname].atomid << std::endl;
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

double get_q(const amoeba_parm &atomid, const t_bonded &bonds, const std::vector<amoeba_parm> &atomids){
    /* If there's only one match, just get the charge */
    int nmatches = atomid.qs.size();
    if (nmatches == 1){
        return atomid.qs[0].q;
    }
    /* Figure out what it's bonded to */
    for (int i=0;i<(int)atomid.qs.size();i++){
        int ztyp = atomid.qs[i].bonded_id[0];
        int xtyp = atomid.qs[i].bonded_id[1];
        int ytyp = atomid.qs[i].bonded_id[2];
        /* Here we are going to see if the parameter exists for stuff it's bound to */
        for (int j=0;j<bonds.nb;j++){
            int jt = atomids[bonds.bond[j]-1].atomid;
            if ( jt == ztyp ){
                for (int k=0;k<bonds.nb;k++){
                    int kt = atomids[bonds.bond[k]-1].atomid;
                    if ( kt == xtyp && k != j ){
                        if (ytyp == 0){
                            return atomid.qs[i].q;
                        }
                        for (int l=0;l<bonds.nb;l++){
                            int lt = atomids[bonds.bond[l]-1].atomid;
                            if (lt == ytyp && l != j && l != k){
                                return atomid.qs[i].q;
                            }
                        }
                    }
                }
            }
        }
    }
    /* Here we are add 1-3 connections */
    for (int i=0;i<(int)atomid.qs.size();i++){
        int ztyp = atomid.qs[i].bonded_id[0];
        int xtyp = atomid.qs[i].bonded_id[1];
        int ytyp = atomid.qs[i].bonded_id[2];
        for (int j=0;j<bonds.nb;j++){
            int jt = atomids[bonds.bond[j]-1].atomid;
            if ( jt == ztyp ){
                for (int k=0;k<bonds.ib;k++){
                    int kt = atomids[bonds.near[k]-1].atomid;
                    bool kpath = true;
                    if (kt == xtyp and kpath){
                        if (ytyp == 0){
                            return atomid.qs[i].q;
                        }
                        for (int l=0;l<bonds.ib;l++){
                            int lt = atomids[bonds.near[k]-1].atomid;
                            bool lpath = true;
                            if (lt == ytyp && l != k && lpath){
                                return atomid.qs[i].q;
                            }
                        }
                    }
                }
            }
        }
    }
    /* Here we hit last resorts */
    for (int i=0;i<(int)atomid.qs.size();i++){
        int ztyp = atomid.qs[i].bonded_id[0];
        int xtyp = atomid.qs[i].bonded_id[1];
        int ytyp = atomid.qs[i].bonded_id[2];
        for (int j=0;j<bonds.nb;j++){
            int jt = atomids[bonds.bond[j]-1].atomid;
            if ( jt == ztyp ){
                if (xtyp == 0){
                    return atomid.qs[i].q;
                }
            }
        }
        if (ztyp == 0){
            return atomid.qs[i].q;
        }
    }
    return 0;
}

int gmx_gmx2xyz(int argc, char *argv[])
{
    const char      *desc[] = {
        "\tConvert from GROMACS format to Tinker XYZ format",
    };

    gmx_bool        bVerbose = FALSE;
    int             a1=0,a2=0;
    static char     *exclude = NULL;
    const char      *tpr_file, *top_file, *traj_file, *index_file, *out_file = NULL;
    char            tpr_title[256], xyz[256];
    int             i_index, switch_lo=5, switch_hi=7;
    atom_id         *ind_index;
    char            *gn_index;
    const char      *parm_name;
    t_pargs         pa[] = {
        { "-a1", TRUE, etINT,
            {&a1}, "Starting atom for bond vector--ie: CD in CNC"},
        { "-a2", TRUE, etINT,
            {&a2}, "Ending atom for bond vector--ie: NE in CNC"},
        { "-exclude", TRUE, etSTR,
            {&exclude}, "Specific atoms to exclude from the field calculations.  Atoms pass by -a1 and -a2 are added to this group automatically.  Note that the whole selection string will need to be quoted so that your shell will pass it in as a string."},
        { "-lo", FALSE, etINT,
            {&switch_lo}, "R_lo for switching function (Angstroms)"},
        { "-hi", FALSE, etINT,
            {&switch_hi}, "R_hi for switching function (Angstroms)."},
        { "-v", FALSE, etBOOL,
            {&bVerbose}, "Be slightly more verbose"}
    };
    t_filenm        fnm[] = {
        {efTPS, NULL, NULL, ffREAD},
        {efTRX, NULL, NULL, ffREAD},
        { efTOP, NULL, NULL, ffREAD },
        { efRND, "-a", "amoeba.prm", ffREAD },
        { efXVG, "-o", "chgpen.xvg", ffWRITE },
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
    FILE            *out = NULL;
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
    top_file    = ftp2fn(efTOP, NFILE, fnm);
    init_top(&top);
    parm_name   = opt2fn("-a", NFILE, fnm);
    out_file    = opt2fn("-o", NFILE, fnm);
    
    /* Open inputs */
    read_tps_conf(tpr_file, buffer, &top, &ePBC,
                  &xtop, NULL, box, TRUE);
    sfree(xtop);
    atoms = &top.atoms;
    
    /* Open xvg output */
    out = xvgropen(out_file,"Charge Penetration Correction for each frame","Frame Number","Field[kbt/eA]", oenv);
    std::vector<std::string> legend;
    legend.push_back("frame");
    legend.push_back("total_field");
    
    /* Make sure switch_lo < switch_hi */
    if (switch_lo >= switch_hi){
        gmx_fatal(FARGS, "\n-lo should be smaller than -hi");
    }
    
    /* Make sure -a1 and -a2 are included and increment by -1 to match internal numbering */
    if ( a1<1 || a2<1 ) {
        gmx_fatal(FARGS, "\nAtom numbers -a1 and -a2 defining the bond vector must be specified\n");
    }
    a1--; a2--;
    
    /* If any, get the atom numbers for the excluded atoms */
    std::vector<int> exclude_ndx;
    if (exclude) {
        int n;
        std::stringstream excluded_line(exclude);
        while ( excluded_line >> n ) {
            exclude_ndx.push_back(n-1);
        }
    }
    exclude_ndx.push_back(a1);
    exclude_ndx.push_back(a2);
    fprintf(stderr,"\nWill calculate and remove field contributions due to atoms:\n\t");
    for (int i=0; i<(int)exclude_ndx.size(); i++) {
        fprintf(stderr,"%i%s ",exclude_ndx[i]+1,*top.atoms.atomname[exclude_ndx[i]]);
        /* Add the excluded atoms to the xvg legend */
        std::stringstream key;
        key << exclude_ndx[i]+1 << *top.atoms.atomname[exclude_ndx[i]] << "_field";
        legend.push_back(key.str());
    }
    fprintf(stderr,"\n");
    
    legend.push_back("F-except");
    
    /* conver the legend to a char array and include it in the xvg file */
    int nlegend_items = legend.size();
    std::vector<std::string> mod_legend;
    char ** flegend = new char * [nlegend_items*3];
    int j = 0;
    for (int i=0;i<nlegend_items*3;i++){
        std::string item = legend[i%nlegend_items];
        if (item == "frame" && i == 0){
            flegend[j] = new char [10*sizeof(item.c_str())];
            strcpy(flegend[j],item.c_str());
            j++;
        }
        else if (item != "frame"){
            item.append("_alpha");
            item.append(std::to_string(i/nlegend_items));
            flegend[j] = new char [10*sizeof(item.c_str())];
            mod_legend.push_back(item);
            strcpy(flegend[j],item.c_str());
            j++;
        }
    }
    int nlegend = mod_legend.size();
    int nexclude = exclude_ndx.size() + 2; // adding 2 for total and f-except
    xvgr_legend(out,j,(const char**)flegend, oenv);
    
    /* Get bonding information */
    t_bonded *bonds;
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
    /* Get 1-3 relationship information from bonded information */
    for (int n=0;n<top.atoms.nr;n++){
        int i = n + 1;
        bonds[i].ib = 0;
        for (int j=0;j<bonds[i].nb;j++){
            int jj = bonds[i].bond[j];
//            std::cout << i << " is bonded to " << jj << std::endl;
            for (int k=0;k<bonds[jj].nb;k++){
                int kk = bonds[jj].bond[k];
                if ( kk != i ){
                    bool cont = true;
//                    std::cout << "\t" << jj << " is bonded to " << kk << std::endl;
                    for (int m=0;m<bonds[kk].nb;m++){
                        int mm = bonds[kk].bond[m];
//                        std::cout << "\t\t" << kk << " is bonded to " << mm << std::endl;
                        if ( mm == i ){
                            cont = false;
                        }
                    }
                    if (cont) {
                        bonds[i].near[bonds[i].ib] = kk;
                        bonds[i].ib++;
//                        std::cout << "\t\t\t" << i << " is 1-3 with " << kk << std::endl;
                    }
                }
            }
        }
    }
    
    /* Get index selection */
    printf("Select group to include in field calculations:\n");
    index_file = ftp2fn(efNDX, NFILE, fnm);
    get_index(atoms, index_file, 1, &i_index, &ind_index, &gn_index);
    if (i_index < 1 or i_index > top.atoms.nr ){
        char error_msg[1052];
        sprintf(error_msg,"Index selection has %i atoms, while the topology file has %i atom entries.\n",i_index,top.atoms.nr);
        gmx_fatal(FARGS,error_msg);
    }
    /* If the index does not start at 1 or if the index skips atoms, this will fail.
       It's easier to just check for that happening than fix it given the current
       implementation, so that's what I'm going to do and just tell the user that it
       failed for this reason
     */
    if (ind_index[0]+1 != 1){
        char error_msg[1052];
        sprintf(error_msg,"Please start the .ndx at atom 1.  Matching to AMOEBA charges will fail otherwise.\nThe .ndx file currently starts with %i.\n",ind_index[0]+1);
        gmx_fatal(FARGS,error_msg);
    }
    for (int i=1;i<i_index;i++){
        if ( ind_index[i] != ind_index[i-1]+1){
            char error_msg[1052];
            sprintf(error_msg,"Please use consecutive atoms in the .ndx file.  Matching to AMOEBA charges will fail otherwise.\nThe problem appears to be index number %i, which is %i, while the previous entry is %i.\n",i,ind_index[i]+1,ind_index[i-1]+1);
            gmx_fatal(FARGS,error_msg);
        }
    }
    
    /* Read AMOEBA parameters */
    std::unordered_map<std::string, std::unordered_map<std::string, amoeba_parm> > atoms_types;
    read_amoeba_parameters(parm_name, atoms_types);
    /* Match AMBER parameters to AMOEBA parameters */
    std::vector<amoeba_parm> atomid_order(i_index);
    std::vector<std::string> atomname_order(i_index);
    for (int i=0;i<i_index;i++){
        atomid_order[i] = get_atomid(ind_index[i],top,atoms_types);
        // Also, since the names are the same every frame, just get those once
        std::stringstream name;
        name << *top.atoms.atomname[ind_index[i]] << "." << *top.atoms.resinfo[top.atoms.atom[ind_index[i]].resind].name;
        atomname_order[i] = name.str().c_str();
    }
    
    /* Read first frame and get the atom charges in order */
    std::vector<double> atomq_order(i_index);
    gmx_bool bHaveFirstFrame = read_first_frame(oenv, &status, traj_file, &fr, flags);
    std::cerr << std::endl;
    if (bHaveFirstFrame) {
        set_trxframe_ePBC(&fr,ePBC);
        for (int n=0;n<i_index;n++){
            int i = ind_index[n];
            atomq_order[n] = get_q(atomid_order[n],bonds[i+1],atomid_order);
//            printf("%i %i %.6f\n",i+1,atomid_order[n].atomid,atomq_order[n]);
            for (int j=0;j<bonds[i+1].nb;j++){
                //printf("%7i(%i) ",bonds[i+1].bond[j],atomid_order[bonds[i+1].bond[j]-1].atomid);
                if (bonds[i+1].bond[j] > ind_index[i_index-1]+1){
                    char error_msg[1052];
                    sprintf(error_msg,"\nPlease use complete molecules.  Matching to AMOEBA charges will fail otherwise.\nAtom %i is connected to atom %i, but the .ndx file stops at atom %i.\n",i+1,bonds[i+1].bond[j],ind_index[i_index-1]+1);
                    gmx_fatal(FARGS,error_msg);
                }
            }
        }
    }
    std::cerr << std::endl;
    /* read file and loop through frames */
    int frameN = 0;
    std::string outfile(out_file);
    do {
        std::vector<float> fieldValues (nlegend);
        float cd[3] = { fr.x[a1][0], fr.x[a1][1],fr.x[a1][2] };
        float ne[3] = { fr.x[a2][0], fr.x[a2][1],fr.x[a2][2] };
        float mid[3];
        rvec bv;
        for (int i=0;i<3;i++){
            mid[i] = (cd[i] + ne[i])/2.;
            bv[i] = 10*(ne[i] - cd[i]);
        }

        std::vector<double> sd_chgpen_field1 (3,0);
        std::vector<double> sd_chgpen_field2 (3,0);
        std::vector<double> sd_chgpen_field3 (3,0);
        std::vector<double> scfield (3,0);

        for (int n=0;n<i_index;n++){
            int i = ind_index[n];
            amoeba_parm parms = atomid_order[n];
            double xyz[3];
            double r,r2,rr,rr2,rr3;
            for (int j=0;j<3;j++){
                xyz[j] = 10 * (mid[j] - fr.x[i][j]);
            }
            r2 = dot(xyz,xyz);
            r = sqrt(r2);
            rr = 1. / r;
            rr2 = rr * rr;
            rr3 = rr * rr2;
            double q_switch = chgpen_switch_function(r,switch_lo,switch_hi);
            double cfield[3];
            double chgpen_field1[3];
            double chgpen_field2[3];
            double chgpen_field3[3];
            
            double d_chgpen_field1[3];
            double d_chgpen_field2[3];
            double d_chgpen_field3[3];

            chgpen_field(xyz, r, atomq_order[n], parms.val.nval, parms.val.alpha1, chgpen_field1);
            chgpen_field(xyz, r, atomq_order[n], parms.val.nval, parms.val.alpha2, chgpen_field2);
            chgpen_field(xyz, r, atomq_order[n], parms.val.nval, parms.val.alpha3, chgpen_field3);

            for (int j=0;j<3;j++){
                cfield[j] = rr3 * atomq_order[n] * xyz[j] * cfac;
                /* Scale with the switching function */
                chgpen_field1[j] = chgpen_field1[j] * q_switch + cfield[j] * (1. - q_switch);
                chgpen_field2[j] = chgpen_field2[j] * q_switch + cfield[j] * (1. - q_switch);
                chgpen_field3[j] = chgpen_field3[j] * q_switch + cfield[j] * (1. - q_switch);
                
                // cfield + x = chgpen, solve for x
                d_chgpen_field1[j] = chgpen_field1[j] - cfield[j];
                d_chgpen_field2[j] = chgpen_field2[j] - cfield[j];
                d_chgpen_field3[j] = chgpen_field3[j] - cfield[j];
                
                sd_chgpen_field1[j] += d_chgpen_field1[j];
                sd_chgpen_field2[j] += d_chgpen_field2[j];
                sd_chgpen_field3[j] += d_chgpen_field3[j];
                scfield[j] += cfield[j];
            }
            /*
            if (r<5){
                std::cout << std::setw(15) << " " << std::setw(15) << project_field(bv,&d_chgpen_field1[0]) << std::setw(15) << project_field(bv,&d_chgpen_field2[0]) << std::setw(15) << project_field(bv,&d_chgpen_field3[0]) << std::endl;
                std::cout << std::setw(15) << project_field(bv,&cfield[0]) << std::setw(15) << project_field(bv,&chgpen_field1[0]) << std::setw(15) << project_field(bv,&chgpen_field2[0]) << std::setw(15) << project_field(bv,&chgpen_field3[0]) << std::endl;
            }
            */
            /* Check to see if this is an excluded atom */
            for (int j=0;j<(int)exclude_ndx.size();j++){
                if (exclude_ndx[j] == i){
                    fieldValues[0*nexclude+j+1] = project_field(bv,&d_chgpen_field1[0]);
                    fieldValues[1*nexclude+j+1] = project_field(bv,&d_chgpen_field2[0]);
                    fieldValues[2*nexclude+j+1] = project_field(bv,&d_chgpen_field3[0]);
                }
            }
        }
        fieldValues[0*nexclude] = project_field(bv,&sd_chgpen_field1[0]);
        fieldValues[1*nexclude-1] = fieldValues[0*nexclude];

        fieldValues[1*nexclude] = project_field(bv,&sd_chgpen_field2[0]);
        fieldValues[2*nexclude-1] = fieldValues[1*nexclude];

        fieldValues[2*nexclude] = project_field(bv,&sd_chgpen_field3[0]);
        fieldValues[3*nexclude-1] = fieldValues[2*nexclude];

        for (int i=1;i<nexclude-1;i++){
            fieldValues[1*nexclude-1] -= fieldValues[0*nexclude+i];
            fieldValues[2*nexclude-1] -= fieldValues[1*nexclude+i];
            fieldValues[3*nexclude-1] -= fieldValues[2*nexclude+i];
        }
        
        fprintf(out,"%10i",frameN);
        for (int i=0;i<(int)mod_legend.size();i++){
//            std::cerr << std::setw(4) << i << std::setw(20) << mod_legend[i] << std::setw(20) << fieldValues[i] << std::endl;
            fprintf(out,"%14.6e",fieldValues[i]);
        }
        fprintf(out,"\n");
        frameN++;
    } while(read_next_frame(oenv, status, &fr));
    atomid_order.clear();
    atomname_order.clear();
    std::cerr << "\nThanks for playing!\n";
    return 0;
}
