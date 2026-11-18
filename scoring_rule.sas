****************************************************************;
******        HP TREE (PROC HPSPLIT) SCORING CODE        ******;
****************************************************************;
 
******              LABELS FOR NEW VARIABLES              ******;
LABEL _Node_ = 'Numero nodo';
LABEL _Leaf_ = 'Numero foglia';
LABEL _WARN_ = 'Avvertimenti';
LABEL P_bio_suitab = 'Previsto: bio_suitab';
 
 _WARN_ = ' ';
 
******      TEMPORARY VARIABLES FOR FORMATTED VALUES      ******;
 
******             ASSIGN OBSERVATION TO NODE             ******;
IF NOT MISSING(asthma_sev_prev) AND ((asthma_sev_prev >= 2.05))
 THEN DO;
  IF NOT MISSING(eosinophils) AND ((eosinophils < 0.2998283020531826))
   THEN DO;
    IF NOT MISSING(base_exac) AND ((base_exac >= 1.903))
     THEN DO;
      _Node_ = 6;
      _Leaf_ = 3;
      P_bio_suitab = 1;
    END;
    ELSE DO;
      _Node_ = 5;
      _Leaf_ = 2;
      P_bio_suitab = 0;
    END;
END;
ELSE DO;
  _Node_ = 4;
  _Leaf_ = 1;
  P_bio_suitab = 1;
END;
END;
ELSE DO;
_Node_ = 1;
_Leaf_ = 0;P_bio_suitab = 0;
END;
****************************************************************;
******     END OF HP TREE (PROC HPSPLIT) SCORING CODE    ******;
****************************************************************;
