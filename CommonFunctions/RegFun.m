function [varargout] = RegFun(Ref, Mov)
Dim = size(Ref, 1); 
if ~ismember(Dim, [2 3] )
    error('Dim is wrong!'); 
end
mm = mean(Ref,2);
ms = mean(Mov,2);
Sshifted = bsxfun(@minus, Mov, ms);  % [Mov(1,:)-ms(1); Mov(2,:)-ms(2); Mov(3,:)-ms(3)];
Mshifted = bsxfun(@minus, Ref, mm); % [Ref(1,:)-mm(1); Ref(2,:)-mm(2); Ref(3,:)-mm(3)];
K = Sshifted*Mshifted';
K = K/size(Ref, 2);
[U A V] = svd(K);
R1 = V*U';
if det(R1)<0
    B = eye(Dim);
    B(Dim,Dim) = det(V*U');
    R1 = V*B*U';
end
t1 = mm - R1*ms;
if nargout == 1
    Ref=eye(Dim+1);
    Ref(1:Dim,1:Dim)=R1;
    Ref(1:Dim,Dim+1)=t1;
    varargout{1} = Ref;
end
if nargout == 2
    varargout{1} = R1;
    varargout{2} = t1;
end
end




