function n = L1norm(p)
   N1 = size(p,1);
   N2 = size(p,2);
   if (N1 ~= N2)
      error('invalid grid spacing')
   end
   h = 1/N1;
   n = sum(sum(abs(p)))*h^2;
