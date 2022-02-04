package casekit.nmr.elucidation.model;

import lombok.*;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

@AllArgsConstructor
@NoArgsConstructor
@Getter
@Setter
@ToString
public class Grouping {

    Map<String, Double> tolerances = new HashMap<>();
    Map<String, Map<Integer, Set<Integer>>> groups;
    Map<String, Map<Integer, Integer>> transformedGroups;
}
