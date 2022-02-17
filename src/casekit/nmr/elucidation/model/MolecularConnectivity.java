package casekit.nmr.elucidation.model;

import casekit.nmr.model.Signal;
import lombok.*;

import java.util.List;
import java.util.Map;
import java.util.Set;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
public class MolecularConnectivity {

    private int index; // e.g. index within PyLSD
    private String atomType;
    private Signal signal;
    private int equivalence;
    private boolean pseudo;
    private List<Integer> protonCounts;
    private List<Integer> hybridizations;
    private List<Integer> hsqc;
    private Map<Integer, Integer[]> hmbc;
    private Map<Integer, Integer[]> cosy;
    private Map<String, Map<Integer, Set<Integer>>> forbiddenNeighbors;
    private Map<String, Map<Integer, Set<Integer>>> setNeighbors;
    private List<Integer> fixedNeighbors;
    private List<Integer> groupMembers;
}
