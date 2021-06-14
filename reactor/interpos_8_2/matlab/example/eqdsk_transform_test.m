set_defaults_matlab

eqdskval=readeqdsk('eqdsksigns.36151t0.5000');
[eqdskval_new]=eqdsk_transform(eqdskval,257,257);

pause

contour(eqdskval_new.rmesh,eqdskval_new.zmesh,eqdskval_new.psi',[1e-5 2e-5],'r--')
contour(eqdskval_new.rmesh,eqdskval_new.zmesh,eqdskval_new.psi',-[1e-5 2e-5],'b')
contour(eqdskval_new.rmesh,eqdskval_new.zmesh,eqdskval_new.psi',-[3e-5 4e-5],'b')
