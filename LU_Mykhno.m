function x = LU_Method(A, B)
  n = length(B);
  for k = 1:n
    L(k,k) = 1;
    
    for j = k:n
      sum = 0;
      for s = 1:k-1
          sum = sum + L(k,s) * U(s,j);
      end
      U(k,j) = A(k,j) - sum;
    end
    
    for i = k+1:n
      sum = 0;  
      for s = 1:k-1
          sum = sum + L(i,s) * U(s,k);
      end  
      L(i,k) =(A(i,k) - sum) / U(k,k);
    end
  end
   
  Z(1) = B(1);

  for i = 2:n
    sum = 0;
    for j = 1:i-1
      sum = sum + L(i,j) * Z(j);
    end
    
    Z(i) = B(i) - sum;
  end
  
  x(n) = Z(n) / U(n,n);
  
  for i = n-1:-1:1
    sum = 0;
    
    for j = i+1:n
      sum = sum + U(i,j) * x(j);
    end
    
    x(i) = (Z(i) - sum) / U(i,i);
  end
 end