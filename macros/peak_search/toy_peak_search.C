#include <TRandom3.h> 
#include <Math/ProbFunc.h>
#include <Math/Functor.h>
#include <Math/Factory.h>
#include <TF1.h> 
#include <TCanvas.h>
#include <TH1D.h> 

#include <cmath>


//this macro is designed to test different Profile Likelihood Ratio (PLR) tests on varieties of TOY (made-up, unrealistic) data.

namespace {

    //snr = the likelihood that any event is 'accidental' or not.
    //x is our 'primary' search variable
    constexpr double x_min = 100.;
    constexpr double x_max = 400.; 
    
    constexpr double pi = 3.1415926536;

    //standerd toy 'event'
    struct event_t { 
        double x; //x-value
        double p_coinc; // probability that this event represents a true 'coincidence'   
    };

    //this will be the background event generator 
    double generate_bg_accidental(TRandom3& rand) {

        static constexpr double x0_bg = 250.; 
        static constexpr double sigma_bg = 45.; 

        double x; 
        do {
            x = rand.Gaus()*sigma_bg + x0_bg;
    
            //now, bias it a bit 
            x += 50.*std::sin( pi*(x-x_min)/((x_max-x_min)) ); 

        } while (x < x_min || x > x_max); 

        return x; 
    }

}

int toy_peak_search(const ULong64_t n_events)
{
    auto hist = new TH1D("h", "background;x;n. events", 200, x_min, x_max);

    TRandom3 rand; 

    for (ULong64_t i=0; i<n_events; i++) hist->Fill( generate_bg_accidental(rand) );

    new TCanvas; 
    hist->SetMaximum(hist->GetMaximum()*1.2);
    hist->SetMinimum(0.);
    hist->Draw("E"); 

    return 0; 
}