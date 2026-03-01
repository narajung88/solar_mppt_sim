%Cost function is minimizing MPPT hardware and wiring cost
%Cost of panels doesn't matter to us since number and cost of panels is
%constant

%M = number of MPPTs
%s = number of panels connected in series with a single string
%Cost = M * params.cost_per_mppt  + M * s * params.cost_per_string

%Define acceptable power loss (in terms of percentage?)

%Find theoretical maximum power
    %MPPT boost * voltage per panel depending on its angle on the car
    %Find current from each panel, multiply by smallest voltage (parallel
    %voltage is limited)

%Find power threshold we have to stay above to be acceptable (P_max *
%(1-power_loss)

%define nonlinear constraint: power must stay above threshold

%Define cost function for matlab to optimize
%Provide starting guess based on simulation

%Can solve cost function for a variety of epsilons (tolerated power loss
%percentages)