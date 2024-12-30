function [R, t, it, obj_err, img_err] = objpose2(P, Qp, initR)



TOL = 1E-5;
EPSILON = 1E-8;

METHOD = 'SVD';

n = size(P,2);

pbar = sum(P,2)/n;
for i = 1:n
  P(:,i) = P(:,i)-pbar;
end

Q(1:3,n) = 0;
for i = 1 : n
  Q(:,i) = [Qp(:,i);1];
end

F(1:3,1:3,1:n) = 0;
V(1:3) = 0;
for i = 1:n
  V = Q(:,i)/Q(3,i);
  F(:,:,i) = (V*V.')/(V.'*V);
end


tFactor = inv(eye(3)-sum(F,3)/n)/n;

it = 0;

Ri = initR;
sum_(1:3,1) = 0;
for i = 1:n
sum_ = sum_ + (F(:,:,i)-eye(3))*Ri*P(:,i);
end
ti = tFactor*sum_;

% calculate error
Qi = xform(P, Ri, ti);
old_err = 0;
vec(1:3,1) = 0;
for i = 1 : n
vec = (eye(3)-F(:,:,i))*Qi(:,i);
old_err = old_err + dot(vec,vec);
end
  

% compute next pose estimate
[Ri, ti, Qi, new_err] = abskernel(P, Qi, F, tFactor, METHOD);
it = it + 1; 

while (abs((old_err-new_err)/old_err) > TOL) && (new_err > EPSILON) 

  old_err = new_err;
  
  % compute the optimal estimate of R
  [Ri, ti, Qi, new_err] = abskernel(P, Qi, F, tFactor, METHOD);
  it = it + 1;

  if it > 20
      break
  end
  
end

R = Ri;
t = ti;
obj_err = sqrt(new_err/n);

if (nargout >= 5) % calculate image-space error
  Qproj = xformproj(P, Ri, ti);
  img_err = 0;
  vec(1:3,1) = 0;
  for i = 1:n
    vec = Qproj(i)-Qp(i);
    img_err = img_err + dot(vec,vec);
  end
  img_err = sqrt(img_err/n);
end

%% correct possible reflection w.r.t the projection center
%if t(3) < 0  %% This is a wrong assumption HERE XXX
% R = -R;     %% Need to be checked in each iteration !!
% t = -t;
%end

% get back to original refernce frame
t = t - Ri*pbar;

% end of OBJPOSE


function t = estimate_t( R,G,F,P,n )

 sum_(1:3,1) = 0;
 for i = 1:n
   sum_ = sum_ + F(:,:,i)*R*P(:,i);
 end
 t = G*sum_;


function [R, t, Qout, err2] = abskernel(P, Q, F, G, method)


n = size(P,2);

for i = 1:n
  Q(:,i) = F(:,:,i)*Q(:,i);
end

% compute P' and Q'
pbar = sum(P,2)/n;
qbar = sum(Q,2)/n;
for i = 1:n
  P(:,i) = P(:,i)-pbar;
  Q(:,i) = Q(:,i)-qbar;
end

if method == 'SVD' % use SVD solution
  % compute M matrix
  M(1:3,1:3) = 0;
  for i = 1:n
    M = M+P(:,i)*Q(:,i).';
  end
  
  % calculate SVD of M
  [U,S,V] = svd(M);
  
  % compute rotation matrix R
  R = V*(U.');

%   disp(['det(R)= ' num2str(det(R)) ]);

  if sign(det(R)) == 1,
    t = estimate_t( R,G,F,P,n );

    if t(3) < 0 ,
%        disp(['t_3 = ' num2str(t(3)) ]);
       %% we need to invert the t 
       R=-[V(:,1:2) -V(:,3)]*U.';
       t = estimate_t( R,G,F,P,n );
     
      % S       V       U       R       t	
	
    end
  else 
    R=[V(:,1:2) -V(:,3)]*U.';
    t = estimate_t( R,G,F,P,n );

    if t(3) < 0 , 
%       disp(['t_3 = ' num2str(t(3)) ]);
      %% we need to invert the t 
      R =- V*(U.');
      t = estimate_t( R,G,F,P,n );
    end
  end

  if det(R) < 0 ,
%     R
%     kl
  end

  if t(3) < 0,
%     t
%     kl
  end

elseif method == 'QTN' % use quaternion solution
  % compute M matrix
  A(1:4,1:4) = 0;
  for i = 1:n
    A = A + qmatQ([1;Q(:,i)]).'*qmatW([1;P(:,i)]);
  end
  
  % Find the largest eigenvalue of A 
  eigs_options.disp = 0;
  [V,D] = eigs(A, eye(size(A)), 1, 'LM', eigs_options);
  
  % compute rotation matrix R from the quaternion that
  % corresponds to the largest egienvalue of A
  %% -> this is wrong -> we need to take the largest -> which is no always the first one
%  kl	
  R = quat2mat(V);

  sum_(1:3,1) = 0;
  for i = 1:n
    sum_ = sum_ + F(:,:,i)*R*P(:,i);
  end
  t = G*sum_;
end


Qout = xform(P, R, t);

% calculate error
err2 = 0;
vec(1:3,1) = 0;
for i = 1 : n
  vec = (eye(3)-F(:,:,i))*Qout(:,i);
  err2 = err2 + dot(vec,vec);
end

% end of ABSKERNEL



