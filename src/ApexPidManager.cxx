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
    fDet_LHRS_prl1->fCell_width_x = 0.084875; 
    fDet_LHRS_prl1->fCell_x0 = 0.5714375;
    fDet_LHRS_prl1->fCell_width_y = 0.500; 
    fDet_LHRS_prl1->fCell_y0 = 0.500 + 0.020;

    // PRL2 (LHRS shower) -------------------------------------------------------------------------- 
    fDet_LHRS_prl2->fZ = 4.77; 
    fDet_LHRS_prl2->fN_rows = 2; 
    fDet_LHRS_prl2->fN_cols = 17; 
    fDet_LHRS_prl2->fCells = fDet_LHRS_prl2->GenerateCells(2, 17); 
    // the LHRS shower has two cell ids swapped 
    fDet_LHRS_prl2->fCells[2].id = 3; 
    fDet_LHRS_prl2->fCells[3].id = 2; 
    //
    fDet_LHRS_prl2->fCell_width_x = 0.084875; 
    fDet_LHRS_prl2->fCell_x0 = 0.5714375;
    fDet_LHRS_prl2->fCell_width_y = 0.500; 
    fDet_LHRS_prl2->fCell_y0 = 0.500 - 0.020;

    // pre-shower (RHRS pre-shower) -------------------------------------------------------------------------- 
    fDet_RHRS_ps->fZ = 3.50; 
    fDet_RHRS_ps->fN_rows = 2; 
    fDet_RHRS_ps->fN_cols = 24; 
    fDet_RHRS_ps->fCells = fDet_RHRS_ps->GenerateCells(2, 24); 
    fDet_RHRS_ps->fCell_width_x = -0.06439; 
    fDet_RHRS_ps->fCell_x0 = -0.77268;
    fDet_RHRS_ps->fCell_width_y = 0.500; 
    fDet_RHRS_ps->fCell_y0 = 0.500 + 0.000;

    // shower (RHRS shower) -------------------------------------------------------------------------- 
    fDet_RHRS_sh->fZ = 3.64; 
    fDet_RHRS_sh->fN_rows = 5; 
    fDet_RHRS_sh->fN_cols = 15; 
    fDet_RHRS_sh->fCells = fDet_RHRS_sh->GenerateCells(5, 15); 
    fDet_RHRS_sh->fCell_width_x = -0.09639; 
    fDet_RHRS_sh->fCell_x0 = -0.686680;
    fDet_RHRS_sh->fCell_width_y = 0.160; 
    fDet_RHRS_sh->fCell_y0 = 0.340;
} 

ClassImp(ApexPidManager); 