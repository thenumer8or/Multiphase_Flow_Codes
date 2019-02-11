MYDIR = ./
FCOMPILER = gfortran -c -O3

OBJS = incomp2d.o fftcos.o fftpackt.o

incomp2d:	$(OBJS)
		gfortran -O3 -o incomp2d $(OBJS)

incomp2d.o:	$(MYDIR)incomp2d.f
		$(FCOMPILER) $(MYDIR)incomp2d.f
		
fftcos.o:	$(MYDIR)fftcos.f
		$(FCOMPILER) $(MYDIR)fftcos.f
		
fftpackt.o:     $(MYDIR)fftpackt.f
		$(FCOMPILER) $(MYDIR)fftpackt.f
