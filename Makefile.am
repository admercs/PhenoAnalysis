SUBDIRS=src config scripts analysis/idl_lib analysis/idl
EXTRA_DIST=data doc analysis/R config_yellowstone.sh config_narsarsuaq.sh config_pilatus.sh config_lema.sh

install-data-local:
	$(MKDIR_P) $(DESTDIR)$(datadir)/sites
	$(INSTALL_DATA) $(srcdir)/data/sites/* $(DESTDIR)$(datadir)/sites/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Example1
	$(INSTALL_DATA) $(srcdir)/data/ini/Example1/* $(DESTDIR)$(datadir)/ini/Example1/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Example2
	$(INSTALL_DATA) $(srcdir)/data/ini/Example2/* $(DESTDIR)$(datadir)/ini/Example2/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Example3
	$(INSTALL_DATA) $(srcdir)/data/ini/Example3/* $(DESTDIR)$(datadir)/ini/Example3/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Global
	$(INSTALL_DATA) $(srcdir)/data/ini/Global/* $(DESTDIR)$(datadir)/ini/Global/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Global-Prediction
	$(INSTALL_DATA) $(srcdir)/data/ini/Global-Prediction/* $(DESTDIR)$(datadir)/ini/Global-Prediction/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Global-Prediction-PFT
	$(INSTALL_DATA) $(srcdir)/data/ini/Global-Prediction-PFT/* $(DESTDIR)$(datadir)/ini/Global-Prediction-PFT/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Assimilation
	$(INSTALL_DATA) $(srcdir)/data/ini/Assimilation/* $(DESTDIR)$(datadir)/ini/Assimilation/
	$(MKDIR_P) $(DESTDIR)$(datadir)/ini/Prediction
	$(INSTALL_DATA) $(srcdir)/data/ini/Prediction/* $(DESTDIR)$(datadir)/ini/Prediction/

uninstall-local:
	rm $(DESTDIR)$(datadir)/sites/*
	rmdir $(DESTDIR)$(datadir)/sites
	rm $(DESTDIR)$(datadir)/ini/Example1/*
	rmdir $(DESTDIR)$(datadir)/ini/Example1
	rm $(DESTDIR)$(datadir)/ini/Example2/*
	rmdir $(DESTDIR)$(datadir)/ini/Example2
	rm $(DESTDIR)$(datadir)/ini/Example3/*
	rmdir $(DESTDIR)$(datadir)/ini/Example3
	rm $(DESTDIR)$(datadir)/ini/Global/*
	rmdir $(DESTDIR)$(datadir)/ini/Global
	rm $(DESTDIR)$(datadir)/ini/Global-Prediction/*
	rmdir $(DESTDIR)$(datadir)/ini/Global-Prediction
	rm $(DESTDIR)$(datadir)/ini/Global-Prediction-PFT/*
	rmdir $(DESTDIR)$(datadir)/ini/Global-Prediction-PFT
	rm $(DESTDIR)$(datadir)/ini/Assimilation/*
	rmdir $(DESTDIR)$(datadir)/ini/Assimilation
	rm $(DESTDIR)$(datadir)/ini/Prediction/*
	rmdir $(DESTDIR)$(datadir)/ini/Prediction
	rmdir $(DESTDIR)$(datadir)/ini
