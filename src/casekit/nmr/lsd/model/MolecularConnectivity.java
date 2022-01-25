package casekit.nmr.lsd.model;

import casekit.nmr.model.Signal;
import lombok.*;

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
    private Set<Integer> protonCounts;
    private Set<Integer> hybridizations;
    private Set<Integer> hsqc;
    private Map<Integer, Integer[]> hmbc;
    private Map<Integer, Integer[]> cosy;
    private Map<String, Map<Integer, Set<Integer>>> forbiddenNeighbors;
    private Map<String, Map<Integer, Set<Integer>>> setNeighbors;
    private Set<Integer> fixedNeighbors;
    private Set<Integer> groupMembers;
}
