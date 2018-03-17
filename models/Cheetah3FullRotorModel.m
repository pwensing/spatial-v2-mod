function [model, graphics] = Cheetah3FullRotorModel()

% This function creates both the rigid body model struct and a graphics
% cell array that specifies the cheetah 3 model

% The 0-configuration for the robot is with legs stright down, cheetah
% pointed along the +x axis of the ICS. 
% The inertial coordinates have +z up (i.e. gravity = [0 0 -9.81 m/s^2])
% The body coordinates have +x forward, +y left, +z up
%

model.gc.point = zeros(3,0);
model.gc.body = zeros(1,0);

Nb = 1;

%% Nominal Paramters of Cheetah
bodyWidth = .226;
bodyHeight = .200;
bodyLength = .600;


motorRad   = .06;
legRad     = 0.03;

lo = 0.044;
lc = 0.13;
l1 = .342; % Length of top link
l2 = .342; % Length of bottom link
motorWidth = 2*lo;
                  
hips{1} = [ bodyLength -bodyWidth 0]'/2;
hips{2} = [ bodyLength  bodyWidth 0]'/2;
hips{3} = [-bodyLength -bodyWidth 0]'/2;
hips{4} = [-bodyLength  bodyWidth 0]'/2;



%% Floating Base Link
model.parent(1) = 0;
model.jtype  = repmat({'  '},13,1);
model.Xtree  = repmat({eye(6)},13,1);
model.I      = repmat({zeros(6)},13,1);

model.parent(1) = 0;
model.jtype{1} = 'Fb';
model.Xtree{1} = eye(6);
model.I{1}     = mcI(16,[0 0 0], boxInertia(32, [bodyLength, bodyWidth, bodyHeight]));
model.Xrotor{1}   = eye(6);
model.Irotor{1}   = zeros(6);
model.gr{1} = 0;


model.gc.body(end+1:end+8) = [Nb, Nb, Nb, Nb, Nb, Nb, Nb, Nb];
model.gc.point(:,end+1:end+8) = [hips{1} hips{2} hips{3} hips{4} hips{1} hips{2} hips{3} hips{4}];
for i = 1:4
    model.gc.point(3,i) = bodyHeight/2;
    model.gc.point(3,i+4) = -0.20; 
end
    
graphics{1}.type = 'box';
graphics{1}.pts  = [-bodyLength -bodyWidth -bodyHeight ; ...
                                   bodyLength  bodyWidth  bodyHeight]/2;

graphics{1}.boundCenter = [0 0 0]';
graphics{1}.boundAxes   = [bodyLength bodyWidth*1.7 bodyHeight*1.7]/2*1.3;





Nb = 1;
side_sign = -1;
%% Legs
for i = 1:4

    %% Ab/Ad (Hip Roll)
    Nb = Nb+1;

    model.parent(Nb) = 1;
    model.jtype{Nb}  = 'Rx';
    model.Xtree{Nb}  = plux(eye(3), [0 0 0]);
    model.Xrotor{Nb} = plux(eye(3), [0 -side_sign*lc 0]');

    model.I{Nb}      = mcI(1.5, [0 side_sign*lo 0], boxInertia(1.5, ones(3,1)*motorWidth));
    model.Irot{Nb}   = mcI(.5, [0 0 0], boxInertia(.5, [.25 1 1]'*motorWidth));

    graphics{Nb}.boundCenter = [0 side_sign*lo 0]';
    graphics{Nb}.boundAxes   = [motorWidth motorWidth motorWidth]/2*1.8;

    graphics{Nb}.boundCenterRot = [0 0 0]';
    graphics{Nb}.boundAxesRot   = [motorWidth motorWidth motorWidth]/2*1.8;

    model.gr{Nb}     = 10.62;


    %% Hip Pitch
    Nb = Nb+1;
    model.parent(Nb) = Nb-1;
    model.jtype{Nb}  = 'Ry';
    model.Xtree{Nb}  = plux(rz(pi), [0 side_sign*lo 0]');
    model.Xrotor{Nb} = plux(rz(pi), [0 side_sign*lo 0]');

    model.I{Nb}      = mcI(.5, [0 0 -l1]/2, boxInertia(.5,[legRad*2 legRad*2 l1]) );
    model.Irot{Nb}   = mcI(.5, [0 0 0]/2, boxInertia(.5, [1 .25 1]'*motorWidth) );

    model.gr{Nb} = 10.62;

    graphics{Nb}.boundCenter = [0 0 -l1]'/2;
    graphics{Nb}.boundAxes   = [legRad*2.2 legRad*2.2 l1*1.2]/2*1.8;

    graphics{Nb}.boundCenterRot = [0 0 0]'/2;
    graphics{Nb}.boundAxesRot   = [motorWidth motorWidth motorWidth]/2*1.8;

    %% Knee Pitch
    Nb = Nb+1;
    model.parent(Nb) = Nb-1;
    model.jtype{Nb}  = 'Ry';

    model.Xtree{Nb}  = plux(eye(3), [0 0 -l1]');
    model.Xrotor{Nb}  = plux(eye(3), [0 0  0 ]');

    model.I{Nb}     = mcI(.5, [0 0 -l2/2], boxInertia(.5,[legRad*2 legRad*2 l2]) );
    model.Irot{Nb}  = mcI(.5, [0 0 0    ], boxInertia(.5,[1 .25 1]'*motorWidth)  );
    model.gr{Nb} = 10.62; 

    model.NB   = Nb;

    graphics{Nb}.boundCenter = [0 0 -l2]'/2;
    graphics{Nb}.boundAxes   = [legRad*2 legRad*2 l2*1.2]/2*1.4;


    graphics{Nb}.boundCenterRot = [0 0 0]';
    graphics{Nb}.boundAxesRot   = [motorWidth motorWidth motorWidth]/2*1.8;
    
    side_sign = -1 * side_sign;
end



mT = 0;
for i = 1:length(model.I)
    mT = mT+model.I{i}(6,6);
end
%mT

model.NB   = Nb;

end

function I = boxInertia(mass, x)
    I = (norm(x)^2*eye(3) - diag(x.^2))*mass/12;
end