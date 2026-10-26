\\ Exact polynomial root counts and independently computed class caps.
cc_rat(x) = type(x)=="t_INT" || type(x)=="t_FRAC";
cc_need(b,s) = if(!b,error(s));
cc_vec(v) = {cc_need(type(v)=="t_VEC","vector required"); for(i=1,#v,cc_need(cc_rat(v[i]),"non-rational entry"));};
cc_pos(p,a,b) = {
  cc_need(cc_rat(a),"rational endpoint required");
  cc_need(subst(p,'t,a)>0,"nonpositive endpoint");
  if(type(b)!="t_INFINITY",cc_need(cc_rat(b),"rational endpoint required");cc_need(subst(p,'t,b)>0,"nonpositive endpoint"));
  if(poldegree(p)>0,cc_need(polsturm(p,[a,b])==0,"real root in positive interval"));
};
cc_B(n,S) = {my(C=1-2/S);cc_need(n*C^2<1,"inadmissible cap");floor(n*(1-C^2)/(1-n*C^2));};
cc_A(n,T,S) = {
  my(U,C,P,a);if(T+S<=4,return(0));
  U=(T+S-4)^2/(4*T*S);C=1-2/S;P=(C+U)/(1-U);a=cc_B(n,S)-1;
  if((n-1)*P^2<1,a=min(a,floor((n-1)*(1-P^2)/(1-(n-1)*P^2))));a;
};
cc_exp_upper(x) = {my(s=1,v=1);for(j=1,512,v=v*x/j;s+=v);1/s;};
cc_atan(a,N) = sum(j=0,N-1,(-1)^j/((2*j+1)*a^(2*j+1)));
cc_verify(D) = {
  my(n,K,tau,E,c,G,d,H,M,J=0,alpha,L0=1,L1,p,q,L2,g,kappa,z,rho,pilo,m,R,oddprod);
  cc_need(type(D)=="t_VEC" && #D==8,"invalid certificate structure");
  n=D[1];K=D[2];tau=D[3];E=D[4];c=D[5];G=D[6];d=D[7];H=D[8];
  cc_need(type(n)=="t_INT" && n>=2,"invalid dimension");
  cc_need(type(K)=="t_INT" && K>=1,"invalid degree");
  cc_need(cc_rat(tau) && tau>0,"invalid scale");cc_vec(E);cc_vec(c);cc_vec(G);
  cc_need(cc_rat(d) && d>0 && cc_rat(H) && H>0,"positive rational claims required");
  M=#E-1;cc_need(M>=2 && #G==M && #c==K+1,"invalid lengths");
  cc_need(E[1]==1 && E[M+1]<4,"invalid edge endpoints");
  for(j=1,M,cc_need(E[j+1]>E[j],"unordered edges");cc_need(G[j]>=0,"negative bin bound");if(E[j]==2,J=j));
  cc_need(J>1 && J<=M,"missing interior squared edge 2");
  cc_need(n*(1-2/E[M+1])^2<1,"inadmissible terminal edge");
  alpha=n/2-1;L1=alpha+1-'t;p=c[1]+c[2]*L1;q=c[1]-c[2]*L1;
  for(k=1,K-1,L2=((2*k+alpha+1-'t)*L1-(k+alpha)*L0)/(k+1);p+=c[k+2]*L2;q+=(-1)^(k+1)*c[k+2]*L2;L0=L1;L1=L2);
  cc_pos(q,0,+oo);cc_pos(-p,tau*E[M+1],+oo);
  for(b=1,M,z=G[b]*(1+('t-tau*E[b])/2+('t-tau*E[b])^2/8)-p;cc_pos(z,tau*E[b],tau*E[b+1]));
  g=vector(M+1,b,if(b<=M,G[b]*cc_exp_upper(tau*E[b]/2),0));
  forstep(b=M-1,J,-1,g[b]=max(g[b],g[b+1]));
  kappa=sum(b=J,M,(g[b]-g[b+1])*cc_B(n,E[b+1]));
  for(j=1,J-1,z=g[j]+sum(b=J,M,(g[b]-g[b+1])*cc_A(n,E[j+1],E[b+1]));kappa=max(kappa,z));
  rho=(subst(p,'t,0)+2*(2^n-1)*kappa)/subst(q,'t,0);cc_need(rho>0,"nonpositive bound");
  pilo=16*cc_atan(5,80)-4*cc_atan(239,79);m=n\2;
  if(n%2==0,
    cc_need(d>=tau^m*rho/(8^m*m!),"density claim too small"),
    oddprod=prod(j=0,m,2*j+1);R=2*tau^m*rho/(4^m*oddprod);cc_need(d^2>=R^2*tau/(8*pilo),"density claim too small")
  );
  cc_need(H^n>=(tau/(2*pilo))^n*rho^2,"Hermite claim too small");
  1;
};
