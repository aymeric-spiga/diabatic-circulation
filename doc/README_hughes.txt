-----Compte rendu des programmes employes dans le cadre du Stage de M2 de recherche portant sur la dynamique atmosphérique de Saturne par Hugues Delattre-----

***** Pour toute questions, contacter hdelattre@hotmail.fr *****


 * Ce fichier README.txt comporte une consigne, un inventaire des objects, une rapide description des dossiers et enfin une rapide description de l'utilisation de chaque code. Il est très FORTEMENT conseillé de lire au minimum la partie consigne pour s'épargner la peine de comprendre par où commencer et lire les codes dans le mauvais ordre de manière indigeste ;).


 * CONSIGNE: ------------------------------------------------------------------
  -Il est conseillé de lire le rapport de stage en premier où au moins en comprendre l'idée principale (notamment la méthode employée). 
  -Il est conseillé en suite de prendre le temps de lire le descriptif des différents programmes et leurs buts. 
  -Pour executer un code -> l'ouvrir avec un éditeur (vi, gedit) et changer les paramètres puis l'exécuter avec python ou ipython nom_du_code.py
  -A la fin de l'execution, si vous appuyez sur ENTREE(retur à la ligne) dans le terminal, toutes les figures se ferment.
  -Chaque code repose sur le même plan, au début on choisit le fichier et les paramètres, puis les constantes planétaires, puis les calculs sont effectués avant d'être affichés dans multiples figures.
  -Les figures sont laissés comme d'origine mais il est possible de les modifier en suivant les commentaires.


 * Ci-joint dans ce dossier :--------------------------------------------------

  archive_photo_articles/...
      ...

  Rapport_de_stage_Hugues.pdf

  README.txt

  simulations_nodyn/
      diagfi_saturn_norings.nc
      diagfi_saturn_rings.nc
      temporaire.nc

  simulations_dynamiques/
      Xhistins_170.nc
      Xhistins_195.nc

  scripthub/
      archives/
          modele_dynamico_solomon.py
          modele_murgatroyd.py
          modele_solomon_diapo.py
          modele_vert_solomon.py  
          theorique_murgatroyd.py        
      coeffbp.txt 
      exemple_script_recup/
          recup_dtetadtphi_model.py
          recup_dtetadtz_model.py 
          recup_modele_wind.py
          recup_Q_model.py
          recup_temp_model.py
          recup_teta_model.py
      modele_compare_rings.py
      modele_dynamico_now.py
      modele_dynamico_w.py
      modele_solomon.py
      theorique_solomon.py



*** Explications :-------------------------------------------------------------

  archive_photo_articles contient (presque) tout les articles employés dans le rapport ainsi que les thèse de Sandrine Guerlet et Melody Sylvestre et la majorité des photos employés dans le rapport (Rapport latex) et dans les présentations.

  simulations_nodyn/ contient les simulations du modèle radiatif convectif de Saturne (Guerlet et al.,2014) sur 5 ans pour temporaire, 10 ans pour diagfi_saturn_norings.nc et diagfi_saturn_rings.nc. La différence des deux dernières simulations étant la présence ou non d'anneaux. diagfi_saturn_rings.nc est la simulation la plus employées.

  simulations_dynamiques/ contient les deux simulations de références utilisées dans le cadre du Stage (selectionnée parmis 20 autres). La simulation Xhistins_195.nc ne contient pas de vent vertical W et fut remplacée par Xhistins_170.nc lors du stage. Xhistins_170.nc est la simulation employée par les codes python dynamico.

  scripthub/ contient l'essentiel des codes employés pour reproduire le travail effectué lors du stage.



 * Liste de l'utilité des codes par ordre de conseil d'utilisation :----------

-> theorique_solomon.py : Code test qui réalise un calcul de la circulation diabatique sur une atmosphère terrestre fictive crée de toute pièce. Le champ de température étant extrêmement simplifié, il est conseillé de parfaire ce code si il doit être utilisé.

-> modele_solomon.py : Code principal du stage : Réalise un calcul de la circulation diabatique depuis les sorties d'un fichier du modèle radiatif/convectif (Guerlet et al.,2014) en appliquant la méthode décrite dans (Solomon et al., 1986).

-> modele_compare_rings.py : Réalise les calculs du code modele_solomon.py sur les simulation du modèle radiatif/convectif avec anneaux et sans anneaux et soustrait les deux.

-> modele_dynamico_w.py : Réalise les calculs des termes de forçage ainsi que des vents moyens résiduels à partir des sorties de DYNAMICO.

-> coeffbp.txt est un fichier .txt utilisé par tout les codes dynamico pour calculer le profil d'altitude.

-> modele_dynamico_now.py : Précédente version de modele_dynamico_w.py ne prenant pas en compte W car les simulations originelle ne l'avait pas.

-> archives/ contient des programmes réalisant des figures précises pour le rapport où la présentation orale. Le répertoire contient également des versions inachevés en prenant les formules de (Murgatroyd et al., 1961), si ces travaux sont repris, il faut changer la méthode d'intégration et vérifier si cela converge bien sur des vents réalistes.

-> exemple_script_recup/ contient des exemples de codes pour récupérer des champs sur une simulation DYNAMICO.
