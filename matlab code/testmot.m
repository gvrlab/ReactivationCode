function N= testmot(rvals,len,rt2)

 
  seqlen=[rt2(:).n];
  
  zhigh=0.1:.02:2.5;
  b0=[-4:.01:4];

  diffpdf=zeros(size(b0));
  nbursts=0;
    
  zval=zeros(size(len));
  zmat=zeros(length(rvals));
    
  [nr nc]=size(rvals);
  for n=seqlen(1):seqlen(end)
      
        
    nn=n-seqlen(1)+1;
    
    
    r0=round(rt2(nn).rvals*100000)/100000;
    z0=(r0-mean(r0))/std(r0);
    
    [h]=histc(z0,b0);
    pdf0=h/sum(h);
    
    
    id0=find(len==n);
    rtmp=(round(rvals(id0)*100000)/100000);
    zval(id0)= (rtmp-mean(r0))/std(r0);
    
    ic = ceil(id0/nr);
    ir = id0-(ic-1)*nr;
    for m=1:length(id0)
      zmat(ir(m),ic(m))=zval(id0(m));
    end
    [dn pdf] = calchist(zval(id0),b0);
    
    nbursts=nbursts+dn;      
    %
    diffpdf = nansum([diffpdf; (pdf-pdf0)*dn]);
    %%
    
     
      
    end
    cdf0=cumsum(pdf0);
    dcdf=cumsum(diffpdf)/nbursts;
    
    N.n= nbursts/2;
    
    for nf=1:length(zhigh)
       
      [dum ith]=min(abs(b0-zhigh(nf)));
      
      
      N.hf(nf).chance=1-cdf0(ith);
      N.hf(nf).f=-dcdf(ith)+N.hf(nf).chance;
      N.hf(nf).pval = binocdf(N.hf(nf).f*N.n, N.n, N.hf(nf).chance,'upper');
      
    end
    %
    nth=min(find(cdf0>.95));
    
    N.flags=sum(zmat>b0(nth));
    %
    
    N.b0=b0;
    N.cdf0=cdf0;
    N.cdf=dcdf+cdf0;
    
    
    
   
    
  end
  
 
  

