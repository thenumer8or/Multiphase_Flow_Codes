MYDIR = ./
FCOMPILER = gfortran -c -O3

OBJS = cylinder.o fftcos.o fftpackt.o

cylinder:	$(OBJS)
		gfortran -O3 -o cylinder $(OBJS)

cylinder.o:	$(MYDIR)cylinder.f
		$(FCOMPILER) $(MYDIR)cylinder.f
		
fftcos.o:	$(MYDIR)fftcos.f
		$(FCOMPILER) $(MYDIR)fftcos.f
		
fftpackt.o:     $(MYDIR)fftpackt.f
		$(FCOMPILER) $(MYDIR)fftpackt.f
