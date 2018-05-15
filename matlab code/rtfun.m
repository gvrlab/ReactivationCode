function rt=rtfun(narr)
  % array of seqeunce lengths
  
  
  % # shuffels
  nrep=100000;
  
  for nn=1:length(narr)
    
    n=narr(nn);
    
    % leftward and rightward template
    s2=1:n;
    s3=randperm(n);
    
    for m=1:nrep
      
      %random replay
      s1=randperm(n);
      
      
      % odd and even rvals are for leftward and rightward replay
      [r]=corrcoef(s1,s2);
      r=round(r*10000)/10000;
      rval(2*m-1)=r(1,2);
      [r]=corrcoef(s1,s3);
      r=round(r*10000)/10000;
      rval(2*m)=r(1,2);
      

    end
    
    
    rt(nn).n=n;
    rt(nn).rvals=rval;
    
  end
  