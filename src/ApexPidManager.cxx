#include <ApexPidManager.h> 
#include <ApexPidDetector.h> 

//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
//_____________________________________________________________________________________________________________________
ApexPidManager::ApexPidManager()
{
    //construct all the 
    fDet_LHRS_prl1 = std::unique_ptr<ApexPidDetector>{new ApexPidDetector}; 
    fDet_LHRS_prl2 = std::unique_ptr<ApexPidDetector>{new ApexPidDetector}; 

    fDet_RHRS_ps = std::unique_ptr<ApexPidDetector>{new ApexPidDetector}; 
    fDet_RHRS_sh = std::unique_ptr<ApexPidDetector>{new ApexPidDetector}; 

    //let's set up each PID detector.

    // PRL1 (LHRS pre-shower) -------------------------------------------------------------------------- 
    fDet_LHRS_prl1->fZ = 4.58; 
    fDet_LHRS_prl1->fN_rows = 2; 
    fDet_LHRS_prl1->fN_cols = 17; 
    fDet_LHRS_prl1->fCells = fDet_LHRS_prl1->GenerateCells(2, 17); 
    fDet_LHRS_prl1->fCell_width_x = 0.150; 
    fDet_LHRS_prl1->fCell_x0 = 1.00;
    fDet_LHRS_prl1->fCell_width_y = -0.500; 
    fDet_LHRS_prl1->fCell_y0 = -0.250;

    // PRL2 (LHRS shower) -------------------------------------------------------------------------- 
    fDet_LHRS_prl2->fZ = 4.77; 
    fDet_LHRS_prl2->fN_rows = 2; 
    fDet_LHRS_prl2->fN_cols = 17; 
    fDet_LHRS_prl2->fCells = fDet_LHRS_prl2->GenerateCells(2, 17); 
    // the LHRS shower has two cell ids swapped 
    fDet_LHRS_prl2->fCells[2].id = 3; 
    fDet_LHRS_prl2->fCells[3].id = 2; 
    //
    fDet_LHRS_prl2->fCell_width_x = 0.145; 
    fDet_LHRS_prl2->fCell_x0 = 1.00;
    fDet_LHRS_prl2->fCell_width_y = -0.500; 
    fDet_LHRS_prl2->fCell_y0 = -0.300;

    // pre-shower (RHRS pre-shower) -------------------------------------------------------------------------- 
    fDet_RHRS_ps->fZ = 3.50; 
    fDet_RHRS_ps->fN_rows = 2; 
    fDet_RHRS_ps->fN_cols = 24; 
    fDet_RHRS_ps->fCells = fDet_RHRS_ps->GenerateCells(2, 24); 
    fDet_RHRS_ps->fCell_width_x = -0.101625; 
    fDet_RHRS_ps->fCell_x0 = -1.23;
    fDet_RHRS_ps->fCell_width_y = +0.500; 
    fDet_RHRS_ps->fCell_y0 = +0.25;

    // shower (RHRS shower) -------------------------------------------------------------------------- 
    fDet_RHRS_sh->fZ = 3.64; 
    fDet_RHRS_sh->fN_rows = 5; 
    fDet_RHRS_sh->fN_cols = 15; 
    fDet_RHRS_sh->fCells = fDet_RHRS_sh->GenerateCells(5, 15); 
    fDet_RHRS_sh->fCell_width_x = -0.153; 
    fDet_RHRS_sh->fCell_x0 = -1.10;
    fDet_RHRS_sh->fCell_width_y = +0.180; 
    fDet_RHRS_sh->fCell_y0 = +0.37;
} 