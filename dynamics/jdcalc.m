function  dotq = jdcalc( jtyp, q, qd )

% jcalc  joint transform and motion subspace matrices.
% [Xj,S]=jcalc(type,q)  returns the joint transform and motion subspace
% matrices for a joint of the given type.  jtyp can be either a string or a
% structure containing a string-valued field called 'code'.  Either way,
% the string contains the joint type code.  For joints that take
% parameters (e.g. helical), jtyp must be a structure, and it must contain
% a field called 'pars', which in turn is a structure containing one or
% more parameters.  (For a helical joint, pars contains a parameter called
% 'pitch'.)  q is the joint's position variable.

if ischar( jtyp )
  code = jtyp;
else
  code = jtyp.code;
end
switch code
  case {'S'}
    dotq=rqd(q,qd);
  case 'fb'
    dotq(1:4)=rqd(q(1:4),qd(1:3) );
    R = rq(q(1:4));
    dotq(5:7)= R * qd(4:6);   
  otherwise
    dotq = qd;
end
