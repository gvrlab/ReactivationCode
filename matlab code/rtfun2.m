function rt=rtfun2(narr)

  
  %random shuffles for motive similarity (no left, right)
  nrep=100000;
  
  for nn=1:length(narr)
    
    n=narr(nn);
    
    s2=1:n;
    
    for m=1:nrep
      s1=randperm(n);
      
      
      [r]=corrcoef(s1,s2);
      r=round(r*1000)/1000;
      rval(m)=r(1,2);
     
      

    end
    
    
    rt(nn).n=n;
    rt(nn).rvals=rval;
    
  end
  