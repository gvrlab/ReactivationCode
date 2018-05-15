function [Nt]=testsession(rval,len, rvalS, lenS, rt)
  
  % rval len: from checktempseq
  % rvalS lenS: from checktempshuffle
  
  % rt: random shuffles for Null-distribution
  %from rtfun 
 
  seqlen=[rt(:).n];
  
  % renaming for convenience
  lenL=len(2:2:end);
  lenR=len(1:2:end-1);
  lenLS=lenS(2:2:end);
  lenRS=lenS(1:2:end-1);
  
  rvalsR=rval(1:2:end-1);
  rvalsL=rval(2:2:end);
  rvalsRS=rvalS(1:2:end-1);
  rvalsLS=rvalS(2:2:end);
  %
  
  zhigh=0.1:.02:2.5;
  b0=[-4:.01:4];
    
  z0R=zeros(size(lenR));
  z0L=zeros(size(lenL));
  zvalR=zeros(size(lenR));
  zvalL=zeros(size(lenL));

  
  %compute all z scores (SSIs)
  for n=seqlen(1):seqlen(end) 
       
      nn=n-seqlen(1)+1;
     
      %z scores from random shuffels
      r0=round(rt(nn).rvals(2:2:end)*100000)/100000;    
      id0R=find(lenR==n);
      id0L=find(lenL==n);
      rtmpR=(round(rvalsRS(id0R)*100000)/100000);
      rtmpL=(round(rvalsLS(id0L)*100000)/100000);
      z0R(id0R)=(rtmpR-mean(r0))/std(r0);
      z0L(id0L)=(rtmpL-mean(r0))/std(r0);
      
      
      %%
      % zscores from
      id0R=find(lenR==n);
      id0L=find(lenL==n);
      rtmpR=(round(rvalsR(id0R)*100000)/100000);
      rtmpL=(round(rvalsL(id0L)*100000)/100000);
      zvalR(id0R)= (rtmpR -mean(r0))/std(r0);
      zvalL(id0L)= (rtmpL -mean(r0))/std(r0);
      
      
            
    end
    
    % compute cdfs
    for nb=1:length(b0)
      z=b0(nb);
      
      % if one direction is above criterion z
      highflagsS=zeros(size(rvalsRS));
      highflagsS(find(z0R>z))=1;
      highflagsS(find(z0L>z))=1;
      cdf0(nb)=1-sum(highflagsS)/length(lenRS);
      
      highflags=zeros(size(rvalsR));
      highflags(find(zvalR>z))=1;
      highflags(find(zvalL>z))=1;
      cdf(nb)=1-sum(highflags)/length(lenR);
     
    end
    
    Nt.b0=b0;
    Nt.cdf0=cdf0;
    Nt.cdf=cdf;
    Nt.n=length(lenR); %=length(lenL); number of bursts
    %%%
    %%

    
    %compute p values for zhigh
    for nh=1:length(zhigh)
      highflags=zeros(size(rvalsR));
      highflags(find(zvalR>zhigh(nh)))=1;
      highflags(find(zvalL>zhigh(nh)))=1;
    
      %
      Nt.hf(nh).flags=highflags;
    

      [dum ith]=min(abs(b0-zhigh(nh)));
      Nt.hf(nh).chance=1-cdf0(ith);  
      Nt.hf(nh).f=1-cdf(ith);
    
      Nt.hf(nh).pval = binocdf(Nt.hf(nh).f*Nt(ns).n, Nt.n, Nt.hf(nh).chance,'upper');     
       
            
      % throw out everything below 20 pop bursts
      if Nt(ns).n<20
        Nt(ns).n=0;Nt(ns).hf(nh).pval=nan;Nt(ns).hf(nh).f=nan;
        Nt(ns).cdf=nan(size(b0));Nt(ns).hf(nh).chance=nan;
      end
      
      
    end
    %%
    
    %% detect significant replays
    nth=min(find(cdf0>.95));
    highflags=zeros(size(rvalsR));
    if length(nth)>0
      highflags(find(zvalR>b0(nth)))=1;
      highflags(find(zvalL>b0(nth)))=1;
    end
    Nt.sigflags=highflags;
    
  end
  
 
  

