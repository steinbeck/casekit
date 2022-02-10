package casekit.nmr.model;

import lombok.*;

@NoArgsConstructor
@AllArgsConstructor
@Getter
@Setter
@ToString
public class PathLength {

    private int min;
    private int max;
    private String source;
}
